function [data]=chkhdr(data,varargin)
%CHKHDR    Check and update header field/values of SAClab records
%
%    Description: CHKHDR(DATA) does a number of consistency checks on
%     SAClab records in DATA.  Timing fields are updated to match the data
%     for both evenly and unevenly sampled records.  Station and event
%     latitude and longitude are checked/fixed.  Distance and azimuth are
%     calculated if the locations are ok and LCALDA is set to true.
%     DEPMIN, DEPMAX, DEPMEN are updated to match the current data.
%     
%     CHKHDR(DATA,FIELD1,VALUE1,...,FIELDN,VALUEN) allows additional header
%     changes to be made before the consistency checks/fixes/updates.
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
%    Data requirements: See Notes above
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
% - dataless support
% - enum checks
%    iftype
%    iztype
%    idep

% check data structure
error(seischk(data,'dep'))

% change header first
data=ch(data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN LOCATION CHECK SECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab header setup (only for getting undef definition)
version=[data.version].';
vers=unique(version);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% number of records
nrecs=numel(data);

% get lcalda (don't mess with those with lcalda ~= 'true')
[lcalda]=glgc(data,'lcalda');
check_loc=strcmpi(lcalda,'true');

% get event/station location info
[evla,evlo,stla,stlo,gcarc,az,baz,dist]=...
    gh(data,'evla','evlo','stla','stlo','gcarc','az','baz','dist');

% version specific
defevla=false(nrecs,1); defevlo=defevla;
defstla=defevla; defstlo=defevla;
for i=1:nver
    % get version group
    v=(version==vers(i));
    
    % which to check
    check=(v & check_loc);
    
    % which lat lon are defined
    defevla(check)=(evla(check)~=h(i).undef.ntype);
    defevlo(check)=(evlo(check)~=h(i).undef.ntype);
    defstla(check)=(stla(check)~=h(i).undef.ntype);
    defstlo(check)=(stlo(check)~=h(i).undef.ntype);
end

% latitude check - needs to be geodetic (90 to -90)
badevla=(evla(defevla)>90 | evla(defevla)<-90);
evla(defevla(badevla))=h(i).undef.ntype; defevla(badevla)=false;
badstla=(stla(defstla)>90 | stla(defstla)<-90);
stla(defstla(badstla))=h(i).undef.ntype; defstla(badstla)=false;

% longitude fix - set for symmetry (-180 to 180)
evlo(defevlo)=mod(evlo(defevlo),360);
stlo(defstlo)=mod(stlo(defstlo),360);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BEGIN DELTA CHECK SECTION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% leven check/fix
leven=glgc(data,'leven');
[delta,odelta,b]=gh(data,'delta','odelta','b');

% get header info
error(lgcchk('leven',leven))
fals=strcmpi(leven,'false');

% make sure .ind field exists (simpler code)
if(~isfield(data,'ind')); data(1).ind=[]; end

% loop over all records, get npts, dep(min,max,men)
len=zeros(nrecs,1); depmax=len; depmin=len; depmen=len; tlen=len; e=len;
for i=1:nrecs
    len(i)=size(data(i).dep,1);
    depmax(i)=max(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmen(i)=mean(data(i).dep(:));
    tlen(i)=size(data(i).ind,1);
    if(fals(i))
        % switch to evenly spaced if no .t data and delta is defined
        if(tlen(i)==0 && delta(i)~=h(vers==version(i)).undef.ntype)
            leven{i}=true; odelta(i)=h(vers==version(i)).undef.ntype;
        % error if .t and .x not the same length
        elseif(tlen(i)~=len(i))
            error('SAClab:chkhdr:dataInconsistent',...
                ['Number of independent variable data does not ',...
                'match number of dependent variable data for record %d'],i)
        % update b, e, odelta, delta
        else
            b(i)=data(i).ind(1);
            e(i)=data(i).ind(end);
            odelta(i)=data(i).ind(2)-data(i).ind(1); % p1 => p2
            delta(i)=(data(i).ind(end)-data(i).ind(1))/(len(i)-1); % avg
        end
    else
        % check delta
        if(delta(i)==h(vers==version(i)).undef.ntype)
            error('SAClab:chkhdr:deltaUndefined',...
                'DELTA field undefined for evenly spaced record %d',i)
        end
        % clear .t and update e, odelta
        data(i).ind=[];
        odelta(i)=h(vers==version(i)).undef.ntype;
        e(i)=b(i)+(len(i)-1)*delta(i);
    end
end

% update header
data=ch(data,'delta',delta,'odelta',odelta,...
    'npts',len,'b',b,'e',e,'leven',leven,...
    'depmax',depmax,'depmin',depmin,'depmen',depmen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   END DELTA CHECK SECTION    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
