function [data]=chkhdr(data,option,varargin)
%CHKHDR    Check and update header field/values of SAClab records
%
%    Description: CHKHDR(DATA) does a number of consistency checks on
%     SAClab records in DATA.  Timing fields are updated to match the data
%     for both evenly and unevenly sampled records.  Station and event
%     latitude and longitude are checked/fixed.  Distance and azimuth are
%     calculated if the locations are ok and LCALDA is set to true.
%     DEPMIN, DEPMAX, DEPMEN are updated to match the current data.
%
%     CHKHDR(DATA,OPTION) allows choosing which checks are performed.
%     Available options are 'THEWORKS', 'LOCATION', 'TIMING', 'ENUMS'.
%     Multiple options can be selected by passing a cell array of options.
%     LOCATION checks the location fields.  TIMING checks timing fields.
%     ENUMS checks enum fields.  THEWORKS does all the available options.
%     
%     CHKHDR(DATA,OPTION,FIELD1,VALUE1,...,FIELDN,VALUEN) makes header 
%     changes before the consistency checks/fixes/updates.
%
%    Notes:
%     - LEVEN must be TRUE or FALSE (error if not)
%     - Unevenly spaced records with DELTA defined and no timing data are
%       changed to evenly spaced (LEVEN set to true)
%     - Unevenly spaced records must have the same number of independent
%       and dependent points (error if not)
%     - STLA and EVLA must be in the range -90 to 90 (changed to undefined
%       if not).
%     - Distance and azimuth calculations assume WGS-84 ellipsoid and
%       utilize Vicenty's method.  If locations are not defined the
%       distance and azimuth fields (GCARC,AZ,BAZ,DIST) are set to
%       undefined.
%
%    System requirements: Matlab 7
%
%    Header changes: DELTA, ODELTA, NPTS, B, E, LEVEN
%                    DEPMEN, DEPMIN, DEPMAX,
%                    EVLA, EVLO, STLA, STLO,
%                    GCARC, AZ, BAZ, DIST
%
%    Usage: data=chkhdr(data)
%           data=chkhdr(data,'field1',values1,...,'fieldN',valuesN)
%
%    Examples:
%     Update header field and also modify the t0 field with 
%     some arrival data:
%      data=chkhdr(data,'t0',P_arrivaltimes)
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 28, 2008 at 22:50 GMT

% todo:
% - categories
%   - all
%   - version   x
%   - vsdata    x
%   - location  x
%   - spacing   x
%   - timing
%   - enums
% - enum checks
%    iftype
%    iztype
%    idep
% - spacing checks
%    check leven vs iftype

% check option
if(nargin==1 || isempty(option)); option='all'; end
option=lower(option);

% change header before checking
data=ch(data,varargin{:});

% grab header setup (will check data structure too)
[h,vi]=vinfo(data);

% number of records
nrecs=numel(data);


if(~isempty(intersect(option,{'version' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN VERSION CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - to verify version listed in header matches that in struct
    
    if(~isequal([data.version].',gh(data,'nvhdr')))
        error('SAClab:chkhdr:verMismatch',...
            ['Version info is corrupted!\n'...
            'One or more records have inconsistent version info!\n'...
            'Check the output of [data.version].'' and gh(data,''nvhdr'').']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END VERSION CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'location' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN LOCATION CHECK SECTION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - To calculate the fields associated with LCALDA if possible.
    % - To look out for and remove bad lat/lon values
    
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


if(~isempty(intersect(option,{'spacing' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN SPACING CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - make sure leven is true/false
    % - make sure delta is defined and >0
    % - make sure npts is defined and >=0
    
    % get leven and delta and npts
    [delta,npts]=gh(data,'delta','npts');
    
    % check leven is 'true' or 'false'
    error(lgcchk('leven',glgc(data,'leven')));
    
    % check npts>=0
    if(npts<0)
        error('SAClab:chkhdr:negativeNPTS',...
            'NPTS cannot be set negative!');
    end
    
    % check delta>0
    if(delta<=0)
        error('SAClab:chkhdr:negativeDELTA',...
                'DELTA cannot be set negative!');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END SPACING CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'vsdata' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN VSDATA CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - check that data matches header info
    % - update dep*
    
    % make sure 'version' and 'spacing' checks are done
    if(isempty(intersect(option,'all')) ...
            && numel(intersect(option,{'version' 'spacing'}))<2)
        % call yourself
        data=chkhdr(data,{'version' 'spacing'});
    end
    
    % get some header info
    ncmp=gncmp(data);
    [b,e,delta,npts]=gh(data,'b','e','delta','npts');
    iftype=genumdesc(data,'iftype');
    leven=glgc(data,'leven');
    
    % loop through each file
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
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
            % check .dep field exists
            if(~isfield(data,'ind'))
                warning('SAClab:chkhdr:noIND',...
                    ['LEVEN set FALSE requires an independent component\n'...
                    'dataset for record %d but there is none:\n'...
                    'Changing LEVEN to TRUE!'],i);
                leven{i}='true';
            end
            
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
        else
            % clear .ind
            if(isfield(data,'ind'))
                data(i).ind=[];
            end
            
            % update e
            e(i)=b(i)+(npts(i)-1)*delta(i);
            
            % update fields to match
            if(npts(i)>0)
                depmax(i)=max(data(i).dep(:));
                depmin(i)=min(data(i).dep(:));
                depmen(i)=mean(data(i).dep(:));
            end
        end
    end
    
    % update header
    warning('off','SAClab:ch:fieldInvalid')
    data=ch(data,'delta',delta,...
        'npts',npts,'ncmp',ncmp,'b',b,'e',e,'leven',leven,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    warning('on','SAClab:ch:fieldInvalid')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END VSDATA CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'timing' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN TIMING CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % What is the point of this section?
    % - make sure delta, npts, b, e agree
    % - 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END TIMING CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if(~isempty(intersect(option,{'enums' 'all'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN ENUMS CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END ENUMS CHECK SECTION    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


end
