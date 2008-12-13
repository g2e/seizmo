function [data]=checkheader(data,options,varargin)
%CHECKHEADER    Check and fix header values of SEIZMO records
%
%    Description: CHECKHEADER(DATA) does a number of consistency checks on
%     SEIZMO records in DATA.  By default it does all of the checks that
%     can be specified (see next calling option).
%
%     CHECKHEADER(DATA,OPTIONS) allows choosing which checks are performed.
%     The available options are:
%
%      'version' - makes sure the version stored in the header and the
%                  struct match
%      'enums'   - assures filetype ids valid, warns about unknown
%                  reference time ids or dependent component ids, checks
%                  that the field that reference time is synced with is 0
%      'leven'   - checks that sample spacing logical is set correctly
%      'delta'   - checks that sample spacing value is valid (>0)
%      'npts'    - checks that npts is valid (>=0)
%      'timing'  - checks NZ fields, checks spectral records are consistent
%                  in regard to freq vs time, makes sure b,e,npts,delta are
%                  consistent
%      'vsdata'  - checks and updates header to match .dep and .ind data
%      'location'- checks locations are valid and gets relative positioning
%                  (if LCALDA is ~FALSE)
%      'all'     - does all of the above
%
%     Note that most of the options build off of one another.  For instance
%     all checks require the version check to have passed.  Multiple
%     options may be given by passing OPTIONS as a cell array of options.
%     
%     CHECKHEADER(DATA,OPTION,FIELD1,VALUE1,...,FIELDN,VALUEN) makes header 
%     changes before the consistency checks and fixes.
%
%    Notes:
%     - CHECKHEADER is LOVROK aware (overwrite ok logical)
%     - use GET_CHECKHEADER_STATE and SET_CHECKHEADER_STATE to check and
%       change if CHECKHEADER is on or off (turning CHECKHEADER off will
%       speed up most functions but will allow inconsistencies to cause
%       touble) 
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA, NPTS, NCMP, B, E, LEVEN,
%                    NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%                    DEPMEN, DEPMIN, DEPMAX, EVLA, EVLO, STLA, STLO,
%                    GCARC, AZ, BAZ, DIST
%
%    Usage: data=checkheader(data)
%           data=checkheader(data,options)
%           data=checkheader(data,options,'field1',values1,...)
%
%    Examples:
%     Check header fields after modifying
%     the t0 field with some arrival data:
%      data=checkheader(data,[],'t0',P_arrivaltimes)
%
%    See also: listheader, changeheader, getheader

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
%        Dec.  4, 2008 - total rewrite, dropped recursion for better
%                        sequencing of checks, significant variable sharing
%                        to reduce repeated calls, more checks, support for
%                        turning this function on/off through SEIZMO global
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  4, 2008 at 05:30 GMT

% todo:
% - excluding certain options
%   - 'exclude_option'
%   - global list
% - check_cmp
%   - inc 0-180
%   - az  0-180
%   - E => az=90, inc=90
%   - N => az=0, inc=90
%   - Z => az=0, inc=0 or 180
%
% testing
%  - single record test of all features
%  - multi record test of all features

% check input
if(nargin<1)
    error('seizmo:checkheader:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin==1 || isempty(options))
    options='all'; 
end

% check SEIZMO global for quick exit
global SEIZMO
if(isfield(SEIZMO,'CHECKHEADER') && isfield(SEIZMO.CHECKHEADER,'ON')...
        && islogical(SEIZMO.CHECKHEADER.ON)...
        && isscalar(SEIZMO.CHECKHEADER.ON)...
        && ~SEIZMO.CHECKHEADER.ON)
    return;
end

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% change header if needed
if(nargin>2); data=changeheader(data,varargin{:}); end

% grab header setup
[h,vi]=versioninfo(data);

% check option is valid
if(~ischar(options) && ~iscellstr(options))
    error('seizmo:checkheader:badOption','OPTIONS must be strings!');
end
valid_options={'all' 'version' 'enums' 'leven' 'delta' 'npts'...
    'timing' 'vsdata' 'location'};
options=unique(lower(cellstr(options)));
if(~isempty(setdiff(options,valid_options)))
    error('seizmo:checkheader:badOption','Unknown OPTIONS!');
end

% check version is consistent (required)
v=getheader(data,'nvhdr');
if(~isequal([data.version].',v))
    i=find([data.version].'~=v);
    error('seizmo:checkheader:verMismatch',...
        ['Version info is corrupted! NVHDR field for record(s):\n'...
        sprintf('%d ',i)...
        '\nis inconsistent with the version stored in the struct!']);
end
options(strcmp('version',options))=[];
if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end

% don't overwrite headers of records if explicitly unallowed
ok=~strcmp(getlgc(data,'lovrok'),'false');
if(any(~ok))
    warning('seizmo:checkheader:lovrokBlock',...
        ['Records:\n' sprintf('%d ',find(~ok))...
        '\nLOVROK is set to FALSE ==> the header cannot be changed!']);
end

% check location if needed
if(~isempty(intersect(options,{'location' 'all'})))
    data=check_location(data,ok,h,vi);
    options(strcmp('location',options))=[];
    if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end
end

% check enums (required at this point)
[nzhour,nzmin,nzsec,nzmsec]=...
    getheader(data,'nzhour','nzmin','nzsec','nzmsec');
[iftype,iztype,idep]=getenumid(data,'iftype','iztype','idep');
[spectral,xyz]=check_enums(data,iftype,iztype,idep,nzhour,nzmin,nzsec,nzmsec);
options(strcmp('enums',options))=[];
if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end

% check leven (required at this point)
leven=getlgc(data,'leven');
[leven]=check_leven(leven,spectral | xyz);
options(strcmp('leven',options))=[];
if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end

% check delta (required at this point)
delta=getheader(data,'delta');
check_delta(delta);
options(strcmp('delta',options))=[];
if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end

% check npts (required at this point)
npts=getheader(data,'npts');
check_npts(npts);
options(strcmp('npts',options))=[];
if(isempty(options)); set_seizmocheck_state(oldseizmocheckstate); return; end

% check timing (required at this point)
[nzyear,nzjday,b,e,sb,sdelta,nspts]=...
    getheader(data,'nzyear','nzjday','b','e','sb','sdelta','nspts');
[nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=...
    check_absolute_timing(nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec,iftype);
check_spectral_timing(b,e,delta,npts,sdelta,nspts,spectral);
[e,delta]=check_normal_timing(b,e,npts,delta,~leven,~xyz & ~spectral);
options(strcmp('timing',options))=[];
if(isempty(options))
    if(any(ok))
        data(ok)=changeheader(data(ok),'nzyear',nzyear(ok),...
            'nzjday',nzjday(ok),'nzhour',nzhour(ok),'nzmin',nzmin(ok),...
            'nzsec',nzsec(ok),'nzmsec',nzmsec(ok),'e',e(ok),...
            'delta',delta(ok));
    end
    set_seizmocheck_state(oldseizmocheckstate);
    return; 
end

% check vs data (required at this point)
[ncmp]=getncmp(data);
[depmin,depmax,depmen]=getheader(data,'depmin','depmax','depmen');
[data,b,e,delta,npts,ncmp,depmin,depmax,depmen,leven]=...
    check_vsdata(data,b,e,delta,npts,ncmp,depmin,depmax,depmen,leven,spectral,h,vi);
options(strcmp('vsdata',options))=[];
if(isempty(options))
    if(any(ok))
        warning('off','seizmo:changeheader:fieldInvalid');
        data(ok)=changeheader(data(ok),'nzyear',nzyear(ok),...
            'nzjday',nzjday(ok),'nzhour',nzhour(ok),'nzmin',nzmin(ok),...
            'nzsec',nzsec(ok),'nzmsec',nzmsec(ok),'e',e(ok),...
            'delta',delta(ok),'npts',npts(ok),'ncmp',ncmp(ok),...
            'b',b(ok),'leven',leven(ok),'depmax',depmax(ok),...
            'depmin',depmin(ok),'depmen',depmen(ok));
        warning('on','seizmo:changeheader:fieldInvalid');
    end
    set_seizmocheck_state(oldseizmocheckstate);
    return; 
end

% update header
if(any(ok))
    warning('off','seizmo:changeheader:fieldInvalid');
    data(ok)=changeheader(data(ok),'nzyear',nzyear(ok),...
        'nzjday',nzjday(ok),'nzhour',nzhour(ok),'nzmin',nzmin(ok),...
        'nzsec',nzsec(ok),'nzmsec',nzmsec(ok),'e',e(ok),...
        'delta',delta(ok),'npts',npts(ok),'ncmp',ncmp(ok),...
        'b',b(ok),'leven',leven(ok),'depmax',depmax(ok),...
        'depmin',depmin(ok),'depmen',depmen(ok));
    warning('on','seizmo:changeheader:fieldInvalid');
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

function [spectral,xyz]=check_enums(data,iftype,iztype,idep,hour,min,sec,msec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BEGIN ENUM CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create valid enum id lists
validftype={'itime' 'irlim' 'iamph' 'ixy' 'ixyz'};
validztype={'iunkn' 'ib' 'iday' 'io' 'ia' 'it0' 'it1' 'it2' 'it3' 'it4'...
    'it5' 'it6' 'it7' 'it8' 'it9'};
validdep={'iunkn' 'idisp' 'ivel' 'iacc' 'ivolts'};

% check data type
if(any(cellfun(@(x)~isempty(setdiff(x,validftype)),iftype)))
    i=find(cellfun(@(x)~isempty(setdiff(x,validftype)),iftype));
    error('seizmo:checkheader:badFileType',...
        ['IFTYPE field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})]);
end

% check dependent component type
if(any(cellfun(@(x)~isempty(setdiff(x,validdep)),idep)))
    i=find(cellfun(@(x)~isempty(setdiff(x,validdep)),idep));
    warning('seizmo:checkheader:badDepType',...
        ['IDEP field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validdep{:})]);
end

% check reference type
badztype=cellfun(@(x)~isempty(setdiff(x,validztype)),iztype);
if(any(badztype))
    i=find(badztype);
    warning('seizmo:checkheader:badRefType',...
        ['IZTYPE field id unknown for record(s):\n' sprintf('%d ',i)...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})]);
end

% get logicals for spectral and xyz datatypes
xyz=strcmp(iftype,'ixyz');
spectral=(strcmp(iftype,'irlim') | strcmp(iftype,'iamph'));

% check day iztype
dayztype=strcmp('iday',iztype);
if(any(dayztype & (hour~=0 | min~=0 | sec~=0 | msec~=0)))
    i=find(dayztype & (hour~=0 | min~=0 | sec~=0 | msec~=0));
    warning('seizmo:checkheader:badRefTime',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNZHOUR, NZMIN, NZSEC, NZMSEC must be 0 for IZTYPE=IDAY!']);
end

% check reference field is zero
goodztype=(~badztype & ~strcmp('iunkn',iztype) & ~dayztype);
if(any(spectral)); iztype(spectral & strcmp('ib',iztype))={'isb'}; end
if(isscalar(unique(iztype(goodztype))))
    iz=unique(iztype(goodztype));
    izt=getheader(data(goodztype),iz{1}(2:end));
    if(any(izt))
        i=find(izt);
        warning('seizmo:checkheader:badRefTime',...
            ['%s field not set to 0 for record(s):\n' sprintf('%d ',i)...
            '\nIZTYPE field demands it is!'],upper(iz{1}(2:end)));
    end
else
    for i=find(goodztype).';
        if(getheader(data(i),iztype{i}(2:end))~=0)
            warning('seizmo:checkheader:badRefTime',...
            'Record %d: %s field not set to 0 when IZTYPE demands it!',...
                i,upper(iztype{i}(2:end)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    END ENUM CHECK SECTION    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [tru]=check_leven(leven,truereq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BEGIN LEVEN CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% leven must be true/false!
tru=strcmp(leven,'true');
fals=strcmp(leven,'false');
if(~all(tru | fals))
    i=find(~(tru | fals));
    error('seizmo:checkheader:levenBad',...
        ['LEVEN field must be set TRUE or FALSE for records:\n'...
        sprintf('%d ',i)]);
end

% leven must be true!
if(any(fals & truereq))
    i=find(fals & truereq);
    error('seizmo:checkheader:levenBad',...
        ['LEVEN field must be set TRUE for Spectral/XYZ records:\n'...
        sprintf('%d ',i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    END LEVEN CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function []=check_delta(delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BEGIN DELTA CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta must be positive
% note that this is only required because negative delta breaks some code
if(any(delta<=0))
    i=find(delta<=0);
    error('seizmo:checkheader:negativeDELTA',...
        ['Records:\n' sprintf('%d ',i) '\nDELTA cannot be set <=0!']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    END DELTA CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function []=check_npts(npts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BEGIN NPTS CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check npts>=0
if(any(npts<0))
    i=find(npts<0);
    error('seizmo:checkheader:negativeNPTS',...
        ['Records:\n' sprintf('%d ',i) '\nNPTS cannot be set negative!']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    END NPTS CHECK SECTION    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [year,jday,hour,minute,sec,msec]=...
    check_absolute_timing(year,jday,hour,minute,sec,msec,iftype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN ABS TIMING CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% which are time series records
ts=strcmp(iftype,'itime');

% check for nonintegers
if(any([year jday hour minute sec msec]...
        ~=fix([year jday hour minute sec msec])))
    i=[year jday hour minute sec msec]...
        ~=fix([year jday hour minute sec msec]);
    warning('seizmo:checkheader:badRefTime',...
        ['Records:\n' sprintf('%d ',i) '\nNZ fields must be integers!']);
    year=fix(year);
    jday=fix(jday);
    hour=fix(hour);
    minute=fix(minute);
    sec=fix(sec);
    msec=fix(msec);
end

% check range ok
badjday=(jday<1 | jday>366 | (jday==366 & ~isleapyear(year)));
badhour=(hour<0 | hour>23);
badmin=(minute<0 | minute>59);
badsec=(sec<0 | sec>60);
badmsec=(msec<0 | msec>999);
if(any(ts & (badjday | badhour | badmin | badsec | badmsec)))
    i=find(ts & (badjday | badhour | badmin | badsec | badmsec));
    warning('seizmo:checkheader:badRefTime',...
        ['Time Series Records:\n' sprintf('%d ',i)...
        '\nReference NZJDAY or NZ time fields not in defined range!\n'...
        'Setting out-of-range NZJDAY to 1 & NZ time fields to 0.']);
    jday(ts & badjday)=1;
    hour(ts & badhour)=0;
    minute(ts & badmin)=0;
    sec(ts & badsec)=0;
    msec(ts & badmsec)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END ABS TIMING CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function []=check_spectral_timing(b,e,d,n,sd,sn,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN SPECTRAL TIMING CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quick exit if none
if(~any(s)); return; end

% get fudge factor
global SEIZMO
try
    fudge=SEIZMO.CHECKHEADER.FUDGE;
catch
    fudge=1e-6;
end

% assure b==0
if(any(b(s)~=0))
    i=find(s & b~=0);
    error('seizmo:checkheader:badB',...
        ['Records:\n' sprintf('%d ',i)...
        '\nB field for Spectral records must be 0!']);
end

% check npts=2^?
if(any(n(s)~=2.^ceil(log2(n(s)))))
    i=find(s & (n~=2.^ceil(log2(n))));
    error('seizmo:checkheader:nptsBAD',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNPTS field for Spectral records must be a power of 2!']);
end

% cross check e is nyquist
if(any(abs(e(s)-1./(2*sd(s)))>fudge*e(s)))
    i=find(s & (abs(e-1./(2*sd))>fudge*e(s)));
    error('seizmo:checkheader:badE',...
        ['Records:\n' sprintf('%d ',i)...
        '\nE field for Spectral records gives the Nyquist frequency\n'...
        'and must be consistent with SDELTA (E=1/(2*SDELTA))!']);
end

% cross check nspts with npts
if(any(sn(s)<0) || any(sn(s)>(n(s)/2)))
    i=find(s & (sn<0 | sn>(n/2)));
    error('seizmo:checkheader:badNPTS',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNSPTS field for Spectral records gives the number of\n'...
        'points in the time domain and must be >0 and <=NPTS/2!']);
end

% cross check delta with npts & sdelta
if(any(abs(d(s)-1./(n(s).*sd(s)))>fudge*d(s)))
    i=find(s & (abs(d-1./(n.*sd))>fudge*d));
    error('seizmo:checkheader:badDELTA',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA field for Spectral records gives the frequency\n'...
        'resolution and must be consistent with SDELTA and NPTS!']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END SPECTRAL TIMING CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [e,delta]=check_normal_timing(b,e,npts,delta,fals,tsxy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BEGIN TIMING CHECK SECTION  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get fudge factor
global SEIZMO
try
    fudge=SEIZMO.CHECKHEADER.FUDGE;
catch
    fudge=1e-6;
end

% for uneven files just adjust delta
delta(fals)=(e(fals)-b(fals))./(npts(fals)-1);

% for even files just adjust e (itime and ixy)
check=~fals & tsxy;
e(check)=b(check)+delta(check).*(npts(check)-1);

% check delta is resolvable for the entire record
d=b+delta-b;
if(any(abs(d-delta)/delta>fudge))
    % require timing to be accurate to 1 millionth of a sample interval
    i=find(abs(d-delta)/delta>fudge);
    warning('seizmo:checkheader:timecorrupt',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA is too insignificant to be well resolved!\n'...
        'This can introduce significant numerical error!\n'...
        'It is recommended to adjust the relative timing so\n'...
        'that the sample interval can be resolved.  This may\n'...
        'require splitting the record into smaller segments\n'...
        'or simply just changing the reference time to be\n'...
        'closer to the time of the record.']);
end
d=e+delta-e;
if(any(abs(d-delta)/delta>fudge))
    % require timing to be accurate to 1 millionth of a sample interval
    i=find(abs(d-delta)/delta>fudge);
    warning('seizmo:checkheader:timecorrupt',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA is too insignificant to be well resolved!\n'...
        'This can introduce significant numerical error!\n'...
        'It is recommended to adjust the relative timing so\n'...
        'that the sample interval can be resolved.  This may\n'...
        'require splitting the record into smaller segments\n'...
        'or simply just changing the reference time to be\n'...
        'closer to the time of the record.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   END TIMING CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [data,b,e,delta,npts,ncmp,depmin,depmax,depmen,tru,h,vi]=...
    check_vsdata(data,b,e,delta,npts,ncmp,depmin,depmax,depmen,tru,s,h,vi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BEGIN VSDATA CHECK SECTION  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get fudge factor
global SEIZMO
try
    fudge=SEIZMO.CHECKHEADER.FUDGE;
catch
    fudge=1e-6;
end

% hasdata
ok=[data.hasdata].';
idx=find(ok);

% quick exit
if(isempty(idx)); return; end

% get dep data size
[nrows,ncols]=arrayfun(@(x)size(x.dep),data(ok));

% check npts
if(any(nrows~=npts(ok)))
    i=idx(nrows~=npts(ok));
    warning('seizmo:checkheader:nptsInconsistent',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNPTS does not match data!\n'...
        'Changing NPTS to match data!']);
    npts(ok)=nrows;
end

% fix/check ncol for spectral
ncols(s(ok))=ncols(s(ok))/2;
if(any(ncols(s(ok))~=fix(ncols(s(ok)))))
    i=idx(ncols(s(ok))~=fix(ncols(s(ok))));
    error('seizmo:checkheader:badDATA',...
        ['Records:\n' sprintf('%d ',i)...
        '\nSpectral records must have an even number of columns!']);
end

% check ncmp
if(any(ncols~=ncmp(ok)))
    i=idx(ncols~=ncmp(ok));
    warning('seizmo:checkheader:nptsInconsistent',...
        ['Records:\n' sprintf('%d ',i)...
        '\nNCMP does not match data!\n'...
        'Changing NCMP to match data!']);
    ncmp(ok)=ncols;
end

% make sure header can handle multiple components if needed
if(any(ok & ncmp>1))
    mok=ok & ncmp>1;
    verch=~[h(vi).mulcmp.valid];
    if(any(verch(mok)))
        % change version
        goch=verch & mok;
        warning('seizmo:checkheader:versNotMulCmp',...
            ['Records:\n' sprintf('%d ',find(goch))...
            '\nCannot handle multiple components!\n'...
            'Changing to a multi-component version!']);
        data(goch)=...
            changeheader(data(goch),'nvhdr',[h(vi(goch)).mulcmp.altver]);
        [data(goch).version]=deal(h(vi(goch)).mulcmp.altver);
        [h,vi]=versioninfo(data);
    end
end

% work on uneven
fals=~tru;
uidx=find(fals & ok);
if(~isempty(uidx))
    % require .ind field for uneven
    if(~isfield(data,'ind'))
        warning('seizmo:checkheader:noIND',...
            ['Records:\n' sprintf('%d ',uidx)...
            '\nLEVEN set FALSE requires an independent component\n'...
            'dataset but there is none.  Changing LEVEN to TRUE!']);
        fals(ok)=false;
    else
        % get ind data size
        [nirows,nicols]=arrayfun(@(x)size(x.ind),data(uidx));
        
        % switch empty .ind to even
        noind=~(nirows.*nicols) & npts(uidx);
        if(any(noind))
            warning('seizmo:checkheader:noINDdata',...
                ['Records:\n' sprintf('%d ',uidx(noind))...
                'LEVEN set FALSE requires an independent component\n'...
                'dataset but there is none.  Changing LEVEN to TRUE!']);
            fals(uidx(noind))=false;
            nirows(noind)=[];
            nicols(noind)=[];
            uidx(noind)=[];
        end
        
        % make sure only 1 .ind cmp
        if(any(nicols>1))
            error('seizmo:checkheader:badNumCMP',...
                ['Records:\n' sprintf('%d ',find(nicols>1))...
                'Too many independent components (only 1 allowed)!'])
        end
        
        % make .dep & .ind consistent npts
        if(any(nirows~=npts(uidx)))
            bad=nirows~=npts(uidx);
            i=uidx(bad);
            warning('seizmo:checkheader:nptsInconsistent',...
                ['Records:\n' sprintf('%d ',i)...
                'NPTS inconsistent between IND and DEP data!\n'...
                'Truncating longer dataset to length of shorter!']);
            for j=find(bad)';
                if(nirows(j)>npts(uidx(j)))
                    data(uidx(j)).ind=data(uidx(j)).ind(1:npts(uidx(j)));
                else
                    npts(uidx(j))=nirows(j);
                    data(uidx(j)).dep=data(uidx(j)).dep(1:nirows(j),:);
                end
            end
        end
        
        % update uneven fields
        b(uidx(npts(uidx)==0))=nan;
        e(uidx(npts(uidx)==0))=nan;
        b(uidx(npts(uidx)>0))=...
            arrayfun(@(x)x.ind(1),data(uidx(npts(uidx)>0)));
        e(uidx(npts(uidx)>0))=...
            arrayfun(@(x)x.ind(end),data(uidx(npts(uidx)>0)));
        delta(uidx(npts(uidx)>1))=...
            arrayfun(@(x)diff(x.ind([1 end])),...
            data(uidx(npts(uidx)>1)))./npts(uidx(npts(uidx)>1));
    end
end

% work on even
tru=~fals;
if(~isempty(tru))
    % clear .ind
    if(isfield(data,'ind'))
        for i=find(tru)'
            data(i).ind=[];
        end
    end
    
    % update e
    tts=tru & ~s;
    e(tts)=b(tts)+(npts(tts)-1).*delta(tts);
end

% update dep stats if npts>0
somepts=npts>0 & ok;
depmin(somepts)=arrayfun(@(x)min(x.dep(:)),data(somepts));
depmax(somepts)=arrayfun(@(x)max(x.dep(:)),data(somepts));
depmen(somepts)=arrayfun(@(x)mean(x.dep(:)),data(somepts));

% check delta is resolvable for the entire record
d=b+delta-b;
if(any(abs(d-delta)/delta>fudge))
    % require timing to be accurate to 1 millionth of a sample interval
    i=find(abs(d-delta)/delta>fudge);
    warning('seizmo:checkheader:timecorrupt',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA is too insignificant to be well resolved!\n'...
        'This can introduce significant numerical error!\n'...
        'It is recommended to adjust the relative timing so\n'...
        'that the sample interval can be resolved.  This may\n'...
        'require splitting the record into smaller segments\n'...
        'or simply just changing the reference time to be\n'...
        'closer to the time of the record.']);
end
d=e+delta-e;
if(any(abs(d-delta)/delta>fudge))
    % require timing to be accurate to 1 millionth of a sample interval
    i=find(abs(d-delta)/delta>fudge);
    warning('seizmo:checkheader:timecorrupt',...
        ['Records:\n' sprintf('%d ',i)...
        '\nDELTA is too insignificant to be well resolved!\n'...
        'This can introduce significant numerical error!\n'...
        'It is recommended to adjust the relative timing so\n'...
        'that the sample interval can be resolved.  This may\n'...
        'require splitting the record into smaller segments\n'...
        'or simply just changing the reference time to be\n'...
        'closer to the time of the record.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   END VSDATA CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [data]=check_location(data,ok,h,vi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN LOCATION CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% get event/station location info
[evla,evlo,stla,stlo,gcarc,az,baz,dist]=...
    getheader(data,'evla','evlo','stla','stlo','gcarc','az','baz','dist');

% check for undefined (version specific)
defevla=false(size(ok)); defevlo=defevla;
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

% get lcalda (don't mess with those with lcalda == 'false')
check_loc=~strcmpi(getlgc(data,'lcalda'),'false');
if(any(check_loc))
    % undefine fields to be calculated
    gcarc(check_loc)=nan; az(check_loc)=nan;
    baz(check_loc)=nan; dist(check_loc)=nan;
    
    % delaz calc (only if all lat lon defined)
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
if(any(ok))
    data(ok)=changeheader(data(ok),'evla',evla(ok),'evlo',evlo(ok),...
        'stla',stla(ok),'stlo',stlo(ok),...
        'gcarc',gcarc(ok),'az',az(ok),'baz',baz(ok),'dist',dist(ok));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  END LOCATION CHECK SECTION  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
