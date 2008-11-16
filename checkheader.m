function [data]=chkhdr(data,option,varargin)
%CHKHDR    Check and update header field/values of SEIZMO records
%
%    Description: CHKHDR(DATA) does a number of consistency checks on
%     SEIZMO records in DATA.  By default it does all of the checks
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
%        June 15, 2008 - doc update
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
% - combine variable space
% - subfunction each section
% - vectorize everything

% check input
if(nargin<1)
    error('seizmo:chkhdr:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin==1 || isempty(option))
    option='all'; 
end
option=lower(option);

% check SEIZMO global for quick exit
global SEIZMO
if(isfield(SEIZMO,'CHKHDR') && isfield(SEIZMO.CHKHDR,'SKIP')...
        && islogical(SEIZMO.CHKHDR.SKIP)...
        && isscalar(SEIZMO.CHKHDR.SKIP) && SEIZMO.CHKHDR.SKIP)
    return;
end

% change header if needed before checking
if(nargin>2); data=ch(data,varargin{:}); end

% grab header setup (will check data structure too)
[h,vi]=vinfo(data);

% turn off struct checking
oldseischkstate=get_seischk_state;
set_seischk_state(true);

% number of records
nrecs=numel(data);


if(~isempty(intersect(option,{'version' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN VERSION CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - to verify version listed in header matches that in struct
    
    if(~isequal([data.version].',gh(data,'nvhdr')))
        error('seizmo:chkhdr:verMismatch',...
            ['Version info is corrupted!\n'...
            'One or more records have inconsistent version info!\n'...
            'Check the output of [data.version].'' '...
            'and gh(data,''nvhdr'').']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END VERSION CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'enums' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN ENUMS CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - make sure enum fields iftype, iztype, idep are valid
    
    % make sure 'version' check is done
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'version'}))<2)
        % call yourself
        data=chkhdr(data,{'version'});
    end
    
    % get enum ids
    [iftype,iztype,idep]=genum(data,'iftype','iztype','idep');
    
    % loop through each record
    for i=1:nrecs
        % valid filetypes
        if(~any(strcmpi(iftype(i),{'itime' 'irlim' 'iamph' 'ixy' 'ixyz'})))
            error('seizmo:chkhdr:badFiletype',...
                ['IFTYPE field (filetype id) for record %d must be of '...
                'the following:\n'...
                'itime irlim iamph ixy ixyz'],i);
        end
        
        % valid reference times
        if(~any(strcmpi(iztype(i),{'iunkn' 'ib' 'iday' 'io' 'ia' 'it0'...
            'it1' 'it2' 'it3' 'it4' 'it5' 'it6' 'it7' 'it8' 'it9'})))
            warning('seizmo:chkhdr:badFiletype',...
                ['IZTYPE field (reference time id) for record %d has '...
                'an unknown value: %s !'],i,iztype{i});
        end
        
        % valid dependent component type
        if(~any(strcmpi(idep(i),{'iunkn' 'idisp' 'ivel' 'iacc' 'ivolts'})))
            warning('seizmo:chkhdr:badFiletype',...
                ['IDEP field (dependent component id) for record %d '...
                'has an unknown value: %s !'],i,idep{i});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END ENUMS CHECK SECTION    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'spacing' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN SPACING CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - make sure leven is true/false
    % - make sure delta is defined and >0
    % - make sure npts is defined and >=0
    % - make sure npts is 2^? for irlim, iamph
    % - make sure leven is true for ixyz, irlim, iamph
    
    % make sure enums check is done (includes version)
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'enums'}))<2)
        % call yourself
        data=chkhdr(data,{'enums'});
    end

    % get leven and delta and npts
    [delta,npts]=gh(data,'delta','npts');
    iftype=genum(data,'iftype');
    leven=glgc(data,'leven');
    
    % check leven is 'true' or 'false'
    error(lgcchk('leven',leven));
    
    % check xyz and spectral files have leven 'true'
    check=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph')...
        | strcmpi(iftype,'ixyz');
    if(any(check))
        if(any(strcmpi(leven(check),'false')))
            error('seizmo:chkhdr:levenBAD',...
                'LEVEN field for Spectral or XYZ records must be true!');
        end
    end
    
    % check npts>=0
    if(npts<0)
        error('seizmo:chkhdr:negativeNPTS',...
            'NPTS cannot be set negative!');
    end
    
    % check delta>0
    if(delta<=0)
        error('seizmo:chkhdr:negativeDELTA',...
                'DELTA cannot be set negative!');
    end
    
    % make sure spectral files have npts == 2^?
    check=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    if(any(check))
        % nextpow2 isn't helpful here
        if(any(npts(check)~=2.^ceil(log2(npts(check)))))
            error('seizmo:chkhdr:nptsBAD',...
                'NPTS field for Spectral records must be a power of 2!');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END SPACING CHECK SECTION   %%%
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
        error('seizmo:chkhdr:badTiming',...
            'B field must not be greater than E field!');
    end
    
    % check spectral files are ok
    spectral=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    if(any(spectral))
        if(any(b(spectral)~=0))
            error('seizmo:chkhdr:badB',...
                'B field for Spectral records must be 0!');
        elseif(any(abs(delta(spectral)...
                -1./(npts(spectral).*sdelta(spectral))))>eps)
            error('seizmo:chkhdr:badDELTA',...
                ['DELTA field for Spectral records must be consistent\n'...
                'with SDELTA and NPTS!']);
        elseif(any(abs(e(spectral)...
                -npts(spectral).*delta(spectral)./2)>eps))
            error('seizmo:chkhdr:badE',...
                ['E field for Spectral records must be consistent\n'...
                'with NPTS and DELTA!']);
        elseif(any(nspts(spectral)<0) ||...
                any(nspts(spectral)>npts(spectral)/2))
            error('seizmo:chkhdr:badNPTS',...
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
                warning('seizmo:chkhdr:badIZTYPE',...
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
            error('seizmo:chkhdr:noDEP',...
                'DEP field required for record %d !',i);
        end
        
        % check that data matches npts, ncmp
        [nrows,ncols]=size(data(i).dep);
        if(any(strcmpi(iftype(i),...
                {'Spectral File-Real/Imag' 'Spectral File-Ampl/Phase'})))
            % spectral files have 2 values per component
            ncols=ncols/2;
            if(fix(ncols)~=ncols)
                error('seizmo:chkhdr:badDATA',...
                    ['Data for spectral record %d does not have '...
                    'an even number of columns!'],i);
            end
        end
        if(nrows~=npts(i))
            warning('seizmo:chkhdr:nptsInconsistent',...
                ['NPTS does not match data for record %d !\n'...
                'Changing NPTS to match data! (%d ==> %d)'],...
                i,npts(i),nrows);
            npts(i)=nrows;
        end
        if(ncols~=ncmp(i))
            warning('seizmo:chkhdr:nptsInconsistent',...
                ['NCMP does not match data for record %d !\n'...
                'Changing NCMP to match data! (%d ==> %d)'],...
                i,ncmp(i),ncols);
            ncmp(i)=ncols;
        end
        
        % make sure header can handle multiple components if needed
        if(ncmp(i)>1)
            if(~h(vi(i)).mulcmp.valid)
                % change version
                warning('seizmo:chkhdr:versNotMulCmp',...
                    ['SEIZMO version %d cannot handle multiple '...
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
                warning('seizmo:chkhdr:noIND',...
                    ['LEVEN set FALSE requires an independent component\n'...
                    'dataset for record %d but there is none:\n'...
                    'Changing LEVEN to TRUE!'],i);
                leven{i}='true';
            else
                % get size
                [nipts,nicmp]=size(data(i).ind);
                
                % switch to evenly spaced if no corresponding independent data
                if((nipts==0 || nicmp==0) && npts(i)>0)
                warning('seizmo:chkhdr:noINDdata',...
                        ['LEVEN set FALSE requires an independent component\n'...
                        'dataset for record %d but there is none:\n'...
                        'Changing LEVEN to TRUE!'],i);
                    leven{i}='true';
                end
                
                % allow only 1 ind cmp
                if(nicmp>1)
                    error('seizmo:chkhdr:badNumCMP',...
                        'Too many independent components for record %d!',i);
                end
                
                % assure consistency with dependent data
                if(nipts~=npts(i))
                    warning('seizmo:chkhdr:nptsInconsistent',...
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
    warning('off','seizmo:ch:fieldInvalid');
    data=ch(data,'delta',delta,...
        'npts',npts,'ncmp',ncmp,'b',b,'e',e,'leven',leven,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    warning('on','seizmo:ch:fieldInvalid');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END VSDATA CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'location' 'all'})))
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

% toggle checking back
set_seischk_state(oldseischkstate);

end
