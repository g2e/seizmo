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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 28, 2008 at 22:50 GMT

% todo:
% - categories
%   - theworks
%   - vsdata
%   - location
%   - spacing
%   - timing
%   - enums
%   - version
% - dataless support
% - enum checks
%    iftype
%    iztype
%    idep

% check data structure
error(seischk(data,'dep'))

% check option
if(nargin<2 || isempty(option)); option='theworks'; end
option=lower(option);

% change header first
data=ch(data,varargin{:});

% grab header setup
[h,vi,v,uv,nv]=vinfo(data);

% number of records
nrecs=numel(data);

if(~isempty(intersect(option,{'theworks' 'location'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN LOCATION CHECK SECTION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get lcalda (don't mess with those with lcalda ~= 'true')
    check_loc=strcmpi(glgc(data,'lcalda'),'true');
    
    % get event/station location info
    [evla,evlo,stla,stlo,gcarc,az,baz,dist]=...
        gh(data,'evla','evlo','stla','stlo','gcarc','az','baz','dist');
    
    % version specific
    defevla=false(nrecs,1); defevlo=defevla;
    defstla=defevla; defstlo=defevla;
    for i=1:nv
        % which to check
        check=(vi==i & check_loc);
        
        % which lat lon are defined
        defevla(check)=(evla(check)~=h(i).undef.ntype);
        defevlo(check)=(evlo(check)~=h(i).undef.ntype);
        defstla(check)=(stla(check)~=h(i).undef.ntype);
        defstlo(check)=(stlo(check)~=h(i).undef.ntype);
    end
    
    % latitude check - needs to be geodetic (90 to -90)
    badevla=(evla(defevla)>90 | evla(defevla)<-90);
    evla(defevla(badevla))=nan; defevla(badevla)=false;
    badstla=(stla(defstla)>90 | stla(defstla)<-90);
    stla(defstla(badstla))=nan; defstla(badstla)=false;
    
    % longitude fix - set to positive (0 to 360)
    evlo(defevlo)=mod(evlo(defevlo),360);
    stlo(defstlo)=mod(stlo(defstlo),360);
    
    % longitude fix - set for symmetry (-180 to 180)
    evlo(defevlo(evlo(defevlo)>180))=evlo(defevlo(evlo(defevlo)>180))-360;
    stlo(defstlo(stlo(defstlo)>180))=stlo(defstlo(stlo(defstlo)>180))-360;
    
    % delaz calc (only if all lat lon defined with lcalda 'true')
    def=(defevla & defevlo & defstla & defstlo);
    if(any(def))
        [gcarc(def),az(def),baz(def),dist(def)]=...
            delaz(evla(def),evlo(def),stla(def),stlo(def));
    end
    
    % update header
    data=ch(data,'evla',evla,'evlo',evlo,'stla',stla,'stlo',stlo,...
        'gcarc',gcarc,'az',az,'baz',baz,'dist',dist);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END LOCATION CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(~isempty(intersect(option,{'theworks' 'spacing'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BEGIN SPACING CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get leven and delta and npts
    leven=glgc(data,'leven');
    [delta,npts]=gh(data,'delta','npts');
    
    % check npts>=0
    if(npts<0)
        error('SAClab:chkhdr:negativeNPTS',...
            'NPTS cannot be set negative!');
    end
    
    % check multi-point records
    mulpts=(npts>1);
    if(any(mulpts))
        error(lgcchk('leven',leven(mulpts)));
        if(delta<=0)
            error('SAClab:chkhdr:negativeDELTA',...
                'DELTA cannot be set negative!');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  END SPACING CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(~isempty(intersect(option,{'theworks' 'vsdata'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN VSDATA CHECK SECTION  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make sure fields exist (simpler code)
    if(~isfield(data,'ind')); data(1).ind=[]; end
    if(~isfield(data,'dep')); data(1).dep=[]; end
    
    % get header info
    leven=glgc(data,'leven');
    uneven=strcmpi(leven,'false');
    [delta,odelta,b,e]=gh(data,'delta','odelta','b','e');
    
    % loop over all records, get npts, dep(min,max,men)
    ilen=nan(nrecs,1); dlen=ilen; ndcmp=ilen; nicmp=ilen;
    depmax=ilen; depmin=ilen; depmen=ilen;
    for i=1:nrecs
        % size up data and get range
        [dlen(i),ndcmp(i)]=size(data(i).dep);
        [ilen(i),nicmp(i)]=size(data(i).ind);
        depmax(i)=max(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmen(i)=mean(data(i).dep(:));
        
        % allow only 1 ind cmp
        if(nicmp(i)>1)
            error('SAClab:chkhdr:badNumCMP',...
                'Too many independent components for record %d!',i);
        end
        
        % check uneven data
        if(uneven(i))
            % switch to evenly spaced if no corresponding independent data
            if(ilen(i)==0 && dlen>0)
                warning('SAClab:chkhdr:levenBad',...
                    ['LEVEN set FALSE without independent '...
                    'dataset for record %d: Changed LEVEN to TRUE!'],i);
                leven{i}=true; odelta(i)=nan;
            % error if .ind and .dep not the same length
            elseif(ilen(i)~=dlen(i))
                error('SAClab:chkhdr:dataInconsistent',...
                    ['Number of independent variable data does not ',...
                    'match the number of dependent variable data for ',...
                    'record %d!'],i)
            % update b, e, odelta, delta
            else
                b(i)=data(i).ind(1);
                e(i)=data(i).ind(end);
                odelta(i)=diff(data(i).ind([1 2])); % pt1 => pt2
                delta(i)=diff(data(i).ind([1 end]))/(dlen(i)-1);
            end
        else
            % check delta
            if(delta(i)>0)
                % clear .ind and update e, odelta
                data(i).ind=[];
                odelta(i)=nan;
                e(i)=b(i)+(dlen(i)-1)*delta(i);
            else
                % clear .ind and update delta, odelta
                warning('SAClab:chkhdr:deltaUnset',...
                    ['DELTA not set for record %d: '...
                    'Using B and E to find DELTA!'],i);
                data(i).ind=[];
                odelta(i)=nan;
                delta(i)=(e(i)-b(i))/(dlen(i)-1);
            end
        end
    end
    
    % update header
    warning('off','SAClab:ch:fieldInvalid')
    data=ch(data,'delta',delta,'odelta',odelta,...
        'npts',dlen,'ncmp',ndcmp,'b',b,'e',e,'leven',leven,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    warning('on','SAClab:ch:fieldInvalid')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END VSDATA CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(~isempty(intersect(option,{'theworks' 'enums'})))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  BEGIN ENUMS CHECK SECTION   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   END ENUMS CHECK SECTION    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end
