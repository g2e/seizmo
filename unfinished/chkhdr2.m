function [data]=chkhdr(data,options,varargin)
%CHKHDR    Check and update header field/values of SAClab records
%
%    Description: CHKHDR(DATA) does a number of consistency checks on
%     SAClab records in DATA.  By default it does all of the checks
%     that can be specified (see next calling option).
%
%     CHKHDR(DATA,OPTIONS) allows choosing which checks are performed.  The
%     available options are:
%
%      'version' - makes sure the version stored in the header and the
%                  struct match
%      'enums'   - assures filetype is valid and warns about a unknown
%                  reference time id or dependent component id
%      'spacing' - checks that delta and npts are valid
%      'timing'  - checks iztype matches with reference field, b<=e, and
%                  updates b,npts,e,delta so they are consistent
%      'vsdata'  - checks and updates header to match data
%      'location'- checks locations are valid and gets relative positioning
%                  if LCALDA is set to TRUE
%      'all'     - does all of the above
%
%     Note that most of the options build off of one another.  For instance
%     all checks require the version check to have passed.  Multiple
%     options may be given by passing OPTIONS as a cell array of options.
%     
%     CHKHDR(DATA,OPTION,FIELD1,VALUE1,...,FIELDN,VALUEN) makes header 
%     changes before the consistency checks/fixes/updates.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA, NPTS, B, E, LEVEN,
%                    DEPMEN, DEPMIN, DEPMAX,
%                    EVLA, EVLO, STLA, STLO,
%                    GCARC, AZ, BAZ, DIST
%
%    Usage: data=chkhdr(data)
%           data=chkhdr(data,options)
%           data=chkhdr(data,options,'field1',values1,...,'fieldN',valuesN)
%
%    Examples:
%     Update header fields and also modify the t0 field with 
%     some arrival data:
%      data=chkhdr(data,[],'t0',P_arrivaltimes)
%
%    See also:  fixdelta, ch, gh

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - renamed to chkhdr, some improvements to
%                        handle unevenly spaced records
%        Feb. 25, 2008 - Couple bugfixes
%        May  12, 2008 - Use new dep* formula
%        June 15, 2008 - Updated documentation
%        June 20, 2008 - Updates header then checks
%        June 22, 2008 - Major revision - supports LCALDA field, better
%                        handling of unevenly sampled data
%        June 28, 2008 - .dep and .ind rather than .x and .t
%        Oct. 15, 2008 - broke into sections, fixed location bugs, now does
%                        location calculations very much like SAC
%        Oct. 26, 2008 - reorganize sections, recursion to allow sections
%                        to build off on another, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 26, 2008 at 21:50 GMT

% todo:
% - this is gonna be used for everything so IT NEEDS TO BE _FAST_!!!
% - combine ch,gh calls
%   - what are we changing/checking?
%       - location - stla,stlo,evla,evlo,gcarc,az,baz,dist,lcalda
%       - version - nvhdr
%       - enums - iftype,iztype,idep
%       - spacing - delta,npts,iftype,leven
%       - timing - b,e,delta,npts,sb,sdelta,nspts,iftype,iztype,leven
%       - vsdata - b,e,delta,npts,leven,iftype,dep*,ncmp
%
%   - how can we arrange these better?
%       - combine variable space
%       - subfunction sections
%       - vectorize everything
%       - break into simpler tasks
%           - by datatype (normal, spectral, uneven, xyz)
%           - by field
%       - no recursion (too many checks)
% - corrective actions vs erroring
%   - maybe an option?
% - excluding certain options
%   - 'exclude_option'
%   - global list

% check input
if(nargin<1)
    error('SAClab:chkhdr:notEnoughInputs','Not enough input arguments.');
elseif(nargin==1 || isempty(options))
    options='all'; 
end

% check SACLAB global for quick exit
global SACLAB
if(isfield(SACLAB,'CHKHDR') && isfield(SACLAB.CHKHDR,'SKIP')...
        && islogical(SACLAB.CHKHDR.SKIP)...
        && isscalar(SACLAB.CHKHDR.SKIP) && SACLAB.CHKHDR.SKIP)
    return;
end

% check data structure
error(seischk(data));

% turn off struct checking
oldseischkstate=get_seischk_state;
set_seischk_state(true);

% change header if needed
if(nargin>2); data=ch(data,varargin{:}); end

% grab header setup
%[h,vi]=vinfo(data);

% unique lowercase options
options=unique(lower(options));

% check option is valid
valid_options={'all' 'version' 'enums' 'spacing'...
    'timing' 'vsdata' 'location'};
if(~ischar(options) || ~iscellstr(options))
    error('SAClab:chkhdr:badOption','OPTIONS must be strings!');
elseif(~isempty(setdiff(options,valid_options)))
    error('SAClab:chkhdr:badOption','Unknown OPTIONS!');
end

% check version is consistent (required)
v=gh(data,'nvhdr');
if(~isequal([data.version].',v))
    i=find([data.version].'~=v);
    error('SAClab:chkhdr:verMismatch',...
        ['Version info is corrupted! NVHDR field for record(s):\n'...
        sprintf('%d ',i)...
        '\nis inconsistent with the version stored in the struct!']);
end
options(strcmp('version',options))=[];
if(isempty(options)); return; end

% check location if needed
if(~isempty(intersect(options,{'location' 'all'})))
    data=check_location(data);
    options(strcmp('location',options))=[];
    if(isempty(options)); return; end
end

% check enums (required at this point)
[nzhour,nzmin,nzsec,nzmsec]=gh(data,'nzhour','nzmin','nzsec','nzmsec');
[iftype,iztype,idep]=genum(data,'iftype','iztype','idep');
[spectral,xyz]=check_enums(iftype,iztype,idep,nzhour,nzmin,nzsec,nzmsec);
options(strcmp('enums',options))=[];
if(isempty(options)); return; end

% check leven (required at this point)
leven=glgc(data,'leven');
[leven]=check_leven(leven,spectral | xyz);
options(strcmp('leven',options))=[];
if(isempty(options)); return; end

% check delta (required at this point)
delta=gh(data,'delta');
check_delta(delta);
options(strcmp('delta',options))=[];
if(isempty(options)); return; end

% check npts (required at this point)
npts=gh(data,'npts');
check_npts(npts);
options(strcmp('npts',options))=[];
if(isempty(options)); return; end

% check timing (required at this point)
[nzyear,nzjday,b,e,sb,sdelta,nspts]=...
    gh(data,'nzyear','nzjday','b','e','sb','sdelta','nspts');
check_absolute_timing(nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec,iftype);
check_spectral_timing(b,e,delta,npts,sb,sdelta,nspts,spectral,iztype);
check_normal_timing(b,e,delta,npts,sb,sdelta,nspts,leven,iftype,iztype);
options(strcmp('timing',options))=[];
if(isempty(options)); return; end

% check vs data (required at this point)
check_vsdata(data,b,e,delta,npts,sb,sdelta,nspts,leven,iftype,iztype);
options(strcmp('vsdata',options))=[];
if(isempty(options)); return; end

% toggle checking back
set_seischk_state(oldseischkstate);

end

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

function [spectral,xyz]=check_enums(iftype,iztype,idep,hour,min,sec,msec)
% check that enum fields iftype, iztype, idep are valid

% create valid enum ids
validftype={'itime' 'irlim' 'iamph' 'ixy' 'ixyz'};
validztype={'iunkn' 'ib' 'iday' 'io' 'ia' 'it0' 'it1' 'it2' 'it3' 'it4'...
    'it5' 'it6' 'it7' 'it8' 'it9'};
validdep={'iunkn' 'idisp' 'ivel' 'iacc' 'ivolts'};

% check data type
if(any(cellfun(@(x)isempty(setdiff(x,validftype)),iftype)))
    i=find(cellfun(@(x)isempty(setdiff(x,validftype)),iftype));
    error('SAClab:chkhdr:badFileType',...
        ['IFTYPE field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validftype)]);
end

% check dependent component type
if(any(cellfun(@(x)isempty(setdiff(x,validdep)),idep)))
    i=find(cellfun(@(x)isempty(setdiff(x,validdep)),idep));
    warning('SAClab:chkhdr:badDepType',...
        ['IDEP field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validdep)]);
end

% check reference type
badztype=cellfun(@(x)isempty(setdiff(x,validztype)),iztype);
if(any(badztype))
    i=find(badztype);
    warning('SAClab:chkhdr:badRefType',...
        ['IZTYPE field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validftype)]);
end

% get logicals for spectral and xyz datatypes
xyz=strcmp(iftype,'ixyz');
spectral=(strcmp(iftype,'irlim') | strcmp(iftype,'iamph'));

% check day iztype
dayztype=strcmp('iday',iztype);
if(any(dayztype & (hour~=0 | min~=0 | sec~=0 | msec~=0)))
    i=find(dayztype & (hour~=0 | min~=0 | sec~=0 | msec~=0));
    warning('SAClab:chkhdr:badRefTime',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNZHOUR, NZMIN, NZSEC, NZMSEC must be 0 for IZTYPE=IDAY!']);
end

% check reference field is zero
goodztype=(~badztype & ~strcmp('iunkn',iztype) & ~dayztype);
if(any(spectral)); iztype(spectral & strcmp('ib',iztype))={'isb'}; end
for i=find(goodztype).';
    if(gh(data(i),iztype{i}(2:end))~=0)
        warning('SAClab:chkhdr:badRefTime',...
            'Record %d: %s field not set to 0 when IZTYPE demands it!',...
            i,iztype{i}(2:end));
    end
end

end

function [tru]=check_leven(leven,truereq)
% check leven is valid

% leven must be true/false!
tru=strcmp(leven,'true');
fals=strcmp(lgc,'false');
if(~all(tru | fals))
    i=find(~(tru | fals));
    error('SAClab:chkhdr:levenBad',...
        ['LEVEN field must be set TRUE or FALSE for records:\n'...
        sprintf('%d ',i)]);
end

% leven must be true!
if(any(fals & truereq))
    i=find(fals & truereq);
    error('SAClab:chkhdr:levenBad',...
        ['LEVEN field must be set TRUE for Spectral/XYZ records:\n'...
        sprintf('%d ',i)]);
end

end

function []=check_delta(delta)
% check delta is valid

% delta must be positive
if(any(delta<=0))
    i=find(delta<=0);
    error('SAClab:chkhdr:negativeDELTA',...
        ['Records:\n' sprintf('%d ',i) '\nDELTA cannot be set negative!']);
end

end

function []=check_npts(npts)
% check npts is valid

% check npts>=0
if(any(npts<0))
    i=find(npts<0);
    error('SAClab:chkhdr:negativeNPTS',...
        ['Records:\n' sprintf('%d ',i) '\nNPTS cannot be set negative!']);
end

end

function [year,jday,hour,minute,sec,msec]=...
    check_absolute_timing(year,jday,hour,minute,sec,msec,iftype)
% check absolute time is valid for time series

% which are time series records
ts=strcmp(iftype,'itime');

% check for nonintegers
if(any([year jday hour minute sec msec]...
        ~=fix([year jday hour minute sec msec])))
    error('SAClab:chkhdr:badRefTime','NZ fields must be integers!');
end

% check range ok
badjday=(jday<1 | jday>366 | (jday==366 & ~isleapyear(year)));
badhour=(hour<0 | hour>23);
badmin=(minute<0 | minute>59);
badsec=(sec<0 | sec>60);
badmsec=(msec<0 | msec>999);
if(any(ts & (badjday | badhour | badmin | badsec | badmsec)))
    i=find(ts & (badjday | badhour | badmin | badsec | badmsec));
    warning('SAClab:chkhdr:badRefTime',...
        ['Time Series Records:\n' sprintf('%d ',i)...
        '\nReference NZJDAY or NZ time fields not in defined range!\n'...
        'Setting out-of-range NZJDAY to 1 & NZ time fields to 0.']);
    jday(ts & badjday)=1;
    hour(ts & badhour)=0;
    minute(ts & badmin)=0;
    sec(ts & badsec)=0;
    msec(ts & badmsec)=0;
end

end

function []=check_spectral_timing(b,e,d,n,sd,sn,s)
% make sure spectral records have frequency/timing consistent

% quick exit if none
if(~any(s)); return; end

% assure b==0
if(any(b(s)~=0))
    i=find(s & b~=0);
    error('SAClab:chkhdr:badB',...
        ['Records:\n' sprintf('%d ',i)...
        '\nB field for Spectral records must be 0!']);
end

% check npts=2^?
if(any(n(s)~=2.^ceil(log2(n(s)))))
    i=find(s & (n~=2.^ceil(log2(n))));
    error('SAClab:chkhdr:nptsBAD',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNPTS field for Spectral records must be a power of 2!']);
end

% cross check e is nyquist
if(any(abs(e(s)-1./(2*sd(s)))>1e-7*e(s)))
    i=find(s & (abs(e-1./(2*sd))>1e-7*e(s)));
    error('SAClab:chkhdr:badE',...
        ['Records:\n' sprintf('%d ',i)...
        '\nE field for Spectral records gives the Nyquist frequency\n'...
        'and must be consistent with SDELTA!']);
end

% cross check nspts with npts
if(any(sn(s)<0) || any(sn(s)>(n(s)/2)))
    i=find(s & (sn<0 | sn>(n/2)));
    error('SAClab:chkhdr:badNPTS',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNSPTS field for Spectral records gives the number of\n'...
        'points in the time domain and must be >0 and less than NPTS/2!']);
end

% cross check delta with npts & sdelta
if(any(abs(d(s)-1./(n(s).*sd(s)))>1e-7*d(s)))
    i=find(s & (abs(d-1./(n.*sd))>1e-7*d));
    error('SAClab:chkhdr:badDELTA',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA field for Spectral records gives the frequency\n'...
        'resolution and must be consistent with SDELTA and NPTS!']);
end

end

function []=check_normal_timing()
% - make sure e>=b
% - make sure delta, npts, b, e agree
% - split by uneven
% - check that e isn't too large for delta
% - only itime?

end

function []=check_vsdata()


end

function [data]=check_location(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN LOCATION CHECK SECTION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - To calculate the fields associated with LCALDA if possible.
    % - To look out for and remove bad lat/lon values
    
    % make sure 'version' check is done
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'version'}))<2)
        % call yourself
        data=chkhdr(data,{'version'});
    end
    
    % get event/station location info
    [evla,evlo,stla,stlo,gcarc,az,baz,dist]=...
        gh(data,'evla','evlo','stla','stlo','gcarc','az','baz','dist');
    
    % check for undefined (version specific)
    defevla=false(nrecs,1); defevlo=defevla;
    defstla=defevla; defstlo=defevla;
    for i=1:numel(h)
        % which to check
        check=(vi==i);
        
        % find lat lon that are defined
        defevla(check)=(evla(check)~=h(i).undef.ntype);
        defevlo(check)=(evlo(check)~=h(i).undef.ntype);
        defstla(check)=(stla(check)~=h(i).undef.ntype);
        defstlo(check)=(stlo(check)~=h(i).undef.ntype);
    end
    
    % latitude check - needs to be 90 to -90
    badevla=(evla>90 | evla<-90);
    evla(defevla & badevla)=nan; 
    defevla(badevla)=false;
    badstla=(stla>90 | stla<-90);
    stla(defstla & badstla)=nan; 
    defstla(badstla)=false;
    
    % longitude fix 
    % - first set to positive (0 to 360)
    % - then set for symmetry (-180 to 180)
    evlo(defevlo)=mod(evlo(defevlo),360);
    stlo(defstlo)=mod(stlo(defstlo),360);
    evlo(defevlo & evlo>180)=evlo(defevlo & evlo>180)-360;
    stlo(defstlo & stlo>180)=stlo(defstlo & stlo>180)-360;
    
    % get lcalda (don't mess with those with lcalda ~= 'true')
    check_loc=strcmpi(glgc(data,'lcalda'),'true');
    if(any(check_loc))
        % undefine fields to be calculated
        gcarc(check_loc)=nan; az(check_loc)=nan;
        baz(check_loc)=nan; dist(check_loc)=nan;
        
        % delaz calc (only if all lat lon defined with lcalda 'true')
        def=(defevla & defevlo & defstla & defstlo & check_loc);
        if(any(def))
            % get geocentric latitude
            geocevla=geodetic2geocentriclat(evla(def));
            geocstla=geodetic2geocentriclat(stla(def));
            
            % get gcarc, az, baz based on sphere (great-circle-arc)
            [gcarc(def),az(def),baz(def)]=...
                sphericalinv(geocevla,evlo(def),geocstla,stlo(def));
            
            % get km dist based on ellipsoid (geodesic)
            dist(def)=vincentyinv(evla(def),evlo(def),stla(def),stlo(def));
        end
    end
    
    % update header
    data=ch(data,'evla',evla,'evlo',evlo,'stla',stla,'stlo',stlo,...
        'gcarc',gcarc,'az',az,'baz',baz,'dist',dist);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END LOCATION CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



if(~isempty(intersect(option,{'timing' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN TIMING CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - make sure e>=b
    % - make sure delta, npts, b, e agree
    % - iztype checks out
    % - watch out for spectral files (b is set 0 rather than -nyq)
    
    % make sure spacing check is done (includes version)
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'spacing'}))<2)
        % call yourself
        data=chkhdr(data,{'spacing'});
    end
    
    % get some header info
    [iftype,iztype]=genum(data,'iftype','iztype');
    leven=glgc(data,'leven');
    [b,e,delta,sb,sdelta,npts,nspts]=...
        gh(data,'b','e','delta','sb','sdelta','npts','nspts');
    
    % e must be greater than b
    if(any(b>e))
        error('SAClab:chkhdr:badTiming',...
            'B field must not be greater than E field!');
    end
    
    % check spectral files are ok
    spectral=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    if(any(spectral))
        if(any(b(spectral)~=0))
            error('SAClab:chkhdr:badB',...
                'B field for Spectral records must be 0!');
        elseif(any(abs(delta(spectral)...
                -1./(npts(spectral).*sdelta(spectral))))>eps)
            error('SAClab:chkhdr:badDELTA',...
                ['DELTA field for Spectral records must be consistent\n'...
                'with SDELTA and NPTS!']);
        elseif(any(abs(e(spectral)...
                -npts(spectral).*delta(spectral)./2)>eps))
            error('SAClab:chkhdr:badE',...
                ['E field for Spectral records must be consistent\n'...
                'with NPTS and DELTA!']);
        elseif(any(nspts(spectral)<0) ||...
                any(nspts(spectral)>npts(spectral)/2))
            error('SAClab:chkhdr:badNPTS',...
                ['NSPTS field for Spectral records must '...
                'be >0 and less than NPTS/2!']);
        end
    end
    
    % check iztype matches up
    check=strcmpi(iztype,'iunkn') | strcmpi(iztype,'ib')...
        | strcmpi(iztype,'iday') | strcmpi(iztype,'io')...
        | strcmpi(iztype,'ia') | strcmpi(iztype,'it0')...
        | strcmpi(iztype,'it1') | strcmpi(iztype,'it2')...
        | strcmpi(iztype,'it3') | strcmpi(iztype,'it4')...
        | strcmpi(iztype,'it5') | strcmpi(iztype,'it6')...
        | strcmpi(iztype,'it7') | strcmpi(iztype,'it8')...
        | strcmpi(iztype,'it9');
    if(any(check))
        for i=find(check).'
            if(gh(data(i),iztype{i}(2:end))~=0)
                warning('SAClab:chkhdr:badIZTYPE',...
                    '%s field not set to 0 when IZTYPE suggests it!',...
                    iztype{i}(2:end));
            end
        end
    end
    
    % make sure b,e,delta,npts matches up
    check=strcmpi(leven,'false');
    delta(check)=(e(check)-b(check))./(npts(check)-1);
    check=~check & ~spectral;
    e(check)=b(check)+delta(check).*(npts(check)-1);
    
    % update header
    data=ch(data,'delta',delta,'e',e);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END TIMING CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'vsdata' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN VSDATA CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - check that data matches header info
    % - update dep*
    
    % make sure timing check is done (includes spacing, version and enums)
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'timing'}))<2)
        % call yourself
        data=chkhdr(data,{'timing'});
    end
    
    % get some header info
    ncmp=gncmp(data);
    [b,e,delta,npts]=gh(data,'b','e','delta','npts');
    [depmen,depmin,depmax]=gh(data,'depmen','depmin','depmax');
    iftype=genumdesc(data,'iftype');
    leven=glgc(data,'leven');
    
    % loop through each file
    for i=1:nrecs
        % skip dataless
        if(~data(i).hasdata); continue; end
        
        % check .dep field exists
        if(~isfield(data,'dep'))
            error('SAClab:chkhdr:noDEP',...
                'DEP field required for record %d !',i);
        end
        
        % check that data matches npts, ncmp
        [nrows,ncols]=size(data(i).dep);
        if(any(strcmpi(iftype(i),...
                {'Spectral File-Real/Imag' 'Spectral File-Ampl/Phase'})))
            % spectral files have 2 values per component
            ncols=ncols/2;
            if(fix(ncols)~=ncols)
                error('SAClab:chkhdr:badDATA',...
                    ['Data for spectral record %d does not have '...
                    'an even number of columns!'],i);
            end
        end
        if(nrows~=npts(i))
            warning('SAClab:chkhdr:nptsInconsistent',...
                ['NPTS does not match data for record %d !\n'...
                'Changing NPTS to match data! (%d ==> %d)'],...
                i,npts(i),nrows);
            npts(i)=nrows;
        end
        if(ncols~=ncmp(i))
            warning('SAClab:chkhdr:nptsInconsistent',...
                ['NCMP does not match data for record %d !\n'...
                'Changing NCMP to match data! (%d ==> %d)'],...
                i,ncmp(i),ncols);
            ncmp(i)=ncols;
        end
        
        % make sure header can handle multiple components if needed
        if(ncmp(i)>1)
            if(~h(vi(i)).mulcmp.valid)
                % change version
                warning('SAClab:chkhdr:versNotMulCmp',...
                    ['SAClab version %d cannot handle multiple '...
                    'components!\nChanging record %d to version %d !'],...
                    data(i).version,i,h(vi(i)).mulcmp.altver);
                data(i)=ch(data(i),'nvhdr',h(vi(i)).mulcmp.altver);
                data(i).version=h(vi(i)).mulcmp.altver;
                [h,vi]=vinfo(data);
            end
        end
        
        % check uneven
        if(strcmpi(leven(i),'false'))
            % check .ind field exists
            if(~isfield(data,'ind'))
                warning('SAClab:chkhdr:noIND',...
                    ['LEVEN set FALSE requires an independent component\n'...
                    'dataset for record %d but there is none:\n'...
                    'Changing LEVEN to TRUE!'],i);
                leven{i}='true';
            else
                % get size
                [nipts,nicmp]=size(data(i).ind);
                
                % switch to evenly spaced if no corresponding independent data
                if((nipts==0 || nicmp==0) && npts(i)>0)
                warning('SAClab:chkhdr:noINDdata',...
                        ['LEVEN set FALSE requires an independent component\n'...
                        'dataset for record %d but there is none:\n'...
                        'Changing LEVEN to TRUE!'],i);
                    leven{i}='true';
                end
                
                % allow only 1 ind cmp
                if(nicmp>1)
                    error('SAClab:chkhdr:badNumCMP',...
                        'Too many independent components for record %d!',i);
                end
                
                % assure consistency with dependent data
                if(nipts~=npts(i))
                    warning('SAClab:chkhdr:nptsInconsistent',...
                        ['NPTS does not match data for record %d !\n'...
                        'Truncating larger to size of smaller!'],i);
                    if(nipts>npts(i))
                        data(i).ind=data(i).ind(1:npts(i));
                    else
                        npts(i)=nipts;
                        data(i).dep=data(i).dep(1:nipts,:);
                    end
                end
                
                % update fields to match
                if(npts(i)==0)
                    b(i)=nan;
                    e(i)=nan;
                else
                    b(i)=data(i).ind(1);
                    e(i)=data(i).ind(end);
                    depmax(i)=max(data(i).dep(:));
                    depmin(i)=min(data(i).dep(:));
                    depmen(i)=mean(data(i).dep(:));
                end
                if(npts(i)>1)
                    delta(i)=diff(data(i).ind([1 end]))/(npts(i)-1);
                end
            end
        end
        if(strcmpi(leven(i),'true'))
            % clear .ind
            if(isfield(data,'ind'))
                data(i).ind=[];
            end
            
            % update e
            e(i)=b(i)+(npts(i)-1)*delta(i);
            
            % fix if spectral record
            if(strcmpi(iftype(i),{'irlim' 'iamph'}))
                e(i)=npts(i)*delta(i)/2;
            end
            
            % update fields to match
            if(npts(i)>0)
                depmax(i)=max(data(i).dep(:));
                depmin(i)=min(data(i).dep(:));
                depmen(i)=mean(data(i).dep(:));
            end
        end
    end
    
    % update header
    warning('off','SAClab:ch:fieldInvalid');
    data=ch(data,'delta',delta,...
        'npts',npts,'ncmp',ncmp,'b',b,'e',e,'leven',leven,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    warning('on','SAClab:ch:fieldInvalid');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END VSDATA CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
