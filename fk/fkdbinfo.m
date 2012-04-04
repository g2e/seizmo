function [maxdb,meddb,mindb]=fkdbinfo(fk,frng,bazrng,srng,esrng,nsrng)
%FKDBINFO    Returns the min/median/max dB for a fk struct
%
%    Usage:    [maxdb,meddb,mindb]=fkdbinfo(fk)
%             [maxdb,meddb,mindb]=fkdbinfo(fk,frng,bazrng,srng,esrng,nsrng)
%
%    Description:
%     [MAXDB,MEDDB,MINDB]=FKDBINFO(FK) returns the decibel limits and
%     median of the beam(s) in the fk struct FK.  This is useful for quick
%     identification for strong plane wave coherency.  The outputs are
%     actually structs with the following format:
%       .db        == decibel value
%       .horzslow  == magnitude of the horizontal slowness (in sec/deg)
%       .backazi   == backazimuth (in degrees)
%       .freq      == frequency (in Hz)
%     Each field is equal in size to the input struct FK.  So MAXDB.db(3),
%     MAXDB.horzslow(3), MAXDB.backazi(3), & MAXDB.freq(3) give info about
%     the max peak in FK(3).  Please note that MEDDB does not return any
%     location info as a median's location is not useful/straightfoward.
%
%     [MAXDB,MEDDB,MINDB]=FKDBINFO(FK,FRNG,BAZRNG,SRNG,ESRNG,NSRNG) returns
%     dB info about a particular subsection of the slowness volume(s) or
%     map(s).  FRNG is the frequency range in Hz.  BAZRNG is the back-
%     azimuth range in degrees.  SRNG is the horizontal slowness magnitude
%     range in sec/deg.  ESRNG & NSRNG are the East & North slowness ranges
%     in sec/deg.  The defaults for all ranges will include all data.
%
%    Notes:
%
%    Examples:
%     % Analyze a 4D fk dataset:
%     s4d=fk4d(data,[],[],50,201,[1/100 1/4]);
%     [maxdb,meddb,mindb]=fkdbinfo(s4d);
%     figure;
%     plot(mindb.db);
%     hold on;
%     plot(maxdb.db);
%     plot(meddb.db);
%     hold off;
%
%    See also: GEOFKDBINFO, FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  26, 2010 - initial version
%        June 30, 2010 - added range arguments
%        July  1, 2010 - range bugfix
%        July  6, 2010 - update for new struct
%        July 12, 2010 - baz range issue bugfix (for sane ranges)
%        July 16, 2010 - output structs with db point info
%        Mar. 29, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 29, 2012 at 14:25 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check fk struct
error(chkfkstruct(fk));

% check ranges
if(nargin<2); frng=[]; end
if(nargin<3); bazrng=[]; end
if(nargin<4); srng=[]; end
if(nargin<5); esrng=[]; end
if(nargin<6); nsrng=[]; end
sf=size(frng);
sbaz=size(bazrng);
ss=size(srng);
ses=size(esrng);
sns=size(nsrng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fkdbinfo:badInput',...
        'FRNG must be a 1x2 array of [FREQLO FREQHI] in Hz!');
end
if(~isempty(bazrng) && (~isreal(bazrng) || numel(sbaz)~=2 || sbaz(1)~=1 ...
        || sbaz(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'BAZRNG must be a 1x2 array of [BAZLO BAZHI] in deg!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:fkdbinfo:badInput',...
        'SRNG must be a 1x2 array of [SLO SHI] in sec/deg!');
end
if(~isempty(esrng) && (~isreal(esrng) || numel(ses)~=2 || ses(1)~=1 ...
        || ses(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'ESRNG must be a 1x2 array of [EAST_SLO EAST_SHI] in sec/deg!');
end
if(~isempty(nsrng) && (~isreal(nsrng) || numel(sns)~=2 || sns(1)~=1 ...
        || sns(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'NSRNG must be a 1x2 array of [NORTH_SLO NORTH_SHI] in sec/deg!');
end

% loop over every element
nfk=numel(fk);
mindb=struct('db',nan(nfk,1,'single'),'horzslow',nan(nfk,1),...
    'backazi',nan(nfk,1),'freq',nan(nfk,1));
meddb=mindb;
maxdb=mindb;
for i=1:nfk
    % get backazimuths and slownesses
    if(fk(i).polar)
        nx=numel(fk(i).x);
        ny=numel(fk(i).y);
        baz=fk(i).x(ones(ny,1),:);
        s=fk(i).y(:,ones(nx,1));
        [es,ns]=slowbaz2kxy(s,baz);
    else % cartesian
        nx=numel(fk(i).x);
        ny=numel(fk(i).y);
        es=fk(i).x(ones(ny,1),:);
        ns=fk(i).y(:,ones(nx,1));
        [s,baz]=kxy2slowbaz(es,ns);
    end
    
    % get extraction indices
    if(fk(i).volume)
        if(isempty(frng))
            fidx=true(1,numel(fk(i).freq));
        else
            fidx=fk(i).freq>=frng(1) & fk(i).freq<=frng(2);
        end
    else
        if(~isempty(frng))
            warning('seizmo:fkdbinfo:badInput',...
                'FRNG does not work for maps!');
        end
        fidx=true;
    end
    if(isempty(bazrng))
        bazidx=true(size(baz));
    else
        bazidx=(baz>=bazrng(1) & baz<=bazrng(2)) ...
            | (baz>=(bazrng(1)-360) & baz<=(bazrng(2)-360)) ...
            | (baz>=(bazrng(1)+360) & baz<=(bazrng(2)+360));
    end
    if(isempty(srng))
        sidx=true(size(s));
    else
        sidx=s>=srng(1) & s<=srng(2);
    end
    if(isempty(esrng))
        esidx=true(size(es));
    else
        esidx=es>=esrng(1) & es<=esrng(2);
    end
    if(isempty(nsrng))
        nsidx=true(size(ns));
    else
        nsidx=ns>=nsrng(1) & ns<=nsrng(2);
    end
    
    % subsection
    idx=repmat(bazidx & sidx & esidx & nsidx,[1 1 numel(fidx)]);
    idx(:,:,~fidx)=false;
    sub=fk(i).beam(idx);
    lind=find(idx);
    [r,c,p]=ind2sub(size(fk(i).beam),lind);
    
    % stats
    % median point cannot be located so it has no point info
    meddb.db(i)=median(sub(:))+fk(i).normdb;
    [mindb.db(i),idx]=min(sub(:));
    mindb.db(i)=mindb.db(i)+fk(i).normdb;
    mindb.horzslow(i)=s(r(idx),c(idx));
    mindb.backazi(i)=baz(r(idx),c(idx));
    mindb.freq(i)=fk(i).freq(p(idx));
    [maxdb.db(i),idx]=max(sub(:));
    maxdb.db(i)=maxdb.db(i)+fk(i).normdb;
    maxdb.horzslow(i)=s(r(idx),c(idx));
    maxdb.backazi(i)=baz(r(idx),c(idx));
    maxdb.freq(i)=fk(i).freq(p(idx));
end

end
