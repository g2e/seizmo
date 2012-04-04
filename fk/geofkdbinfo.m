function [maxdb,meddb,mindb]=geofkdbinfo(fk,frng,srng,latrng,lonrng)
%GEOFKDBINFO    Returns the min/median/max dB for a geofk struct
%
%    Usage:    [maxdb,meddb,mindb]=geofkdbinfo(fk)
%              [maxdb,meddb,mindb]=geofkdbinfo(fk,frng,srng,latrng,lonrng)
%
%    Description:
%     [MAXDB,MEDDB,MINDB]=GEOFKDBINFO(FK) returns the decibel limits and
%     median of the beam(s) in the geofk struct FK.  This is useful for
%     identification of elements with strong spherical wave coherency.  The
%     outputs are actually structs with the following format:
%       .db        == decibel value
%       .horzslow  == magnitude of the horizontal slowness (in sec/deg)
%       .latlon    == [latitude longitude] (in degrees)
%       .freq      == frequency (in Hz)
%     Each field is equal in size to the input struct FK.  So MAXDB.db(3),
%     MAXDB.horzslow(3), MAXDB.latlon(3,:), & MAXDB.freq(3) give info about
%     the max peak in FK(3).  Please note that MEDDB does not return any
%     location info as a median's location is not useful/straightfoward.
%
%     [MAXDB,MEDDB,MINDB]=GEOFKDBINFO(FK,FRNG,SRNG,LATRNG,LONRNG) returns
%     dB info about a particular subsection of the slowness volume(s) or
%     map(s).  FRNG is the frequency range in Hz.  SRNG is the horizontal
%     slowness magnitude range in sec/deg.  LATRNG & LONRNG are the
%     latitude and longitude limits of an area on the Earth in units of
%     degrees.
%
%    Notes:
%
%    Examples:
%
%    See also: FKDBINFO, GEOFKXCVOLUME, GEOFKSUBVOL, GEOFKVOL2MAP

%     Version History:
%        July 12, 2010 - initial version
%        July 18, 2010 - output structs with db point info
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 22:25 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check fk struct
error(chkgeofkstruct(fk));

% check ranges
if(nargin<2); frng=[]; end
if(nargin<3); srng=[]; end
if(nargin<4); latrng=[]; end
if(nargin<5); lonrng=[]; end
sf=size(frng);
ss=size(srng);
sla=size(latrng);
slo=size(lonrng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:geofkdbinfo:badInput',...
        'FRNG must be a 1x2 array of [FREQLO FREQHI] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkdbinfo:badInput',...
        'SRNG must be a 1x2 array of [SLO SHI] in sec/deg!');
end
if(~isempty(latrng) && (~isreal(latrng) || numel(sla)~=2 || sla(1)~=1 ...
        || sla(2)~=2))
    error('seizmo:geofkdbinfo:badInput',...
        'LATRNG must be a 1x2 array of [MIN_LAT MAX_LAT] in degrees!');
end
if(~isempty(lonrng) && (~isreal(lonrng) || numel(slo)~=2 || slo(1)~=1 ...
        || slo(2)~=2))
    error('seizmo:geofkdbinfo:badInput',...
        'LONRNG must be a 1x2 array of [MIN_LON MAX_LON] in degrees!');
end

% loop over every element
nfk=numel(fk);
mindb=struct('db',nan(nfk,1,'single'),'horzslow',nan(nfk,1),...
    'latlon',nan(nfk,2),'freq',nan(nfk,1));
meddb=mindb;
maxdb=mindb;
for i=1:numel(fk)
    % get extraction indices
    if(fk(i).volume(2))
        if(isempty(frng))
            fidx=1:numel(fk(i).freq);
        else
            fidx=find(fk(i).freq>=frng(1) & fk(i).freq<=frng(2));
        end
    else
        if(~isempty(frng))
            warning('seizmo:geofkdbinfo:badInput',...
                'FRNG does not work for frequency maps!');
        end
        fidx=1;
    end
    if(fk(i).volume(1))
        if(isempty(srng))
            sidx=1:numel(fk(i).horzslow);
        else
            sidx=find(fk(i).horzslow>=srng(1) & fk(i).horzslow<=srng(2));
        end
    else
        if(~isempty(srng))
            warning('seizmo:geofkdbinfo:badInput',...
                'SRNG does not work for slowness maps!');
        end
        sidx=1;
    end
    if(isempty(latrng))
        latidx=true(size(fk(i).latlon,1),1);
    else
        latidx=fk(i).latlon>=latrng(1) & fk(i).latlon<=latrng(2);
    end
    if(isempty(lonrng))
        lonidx=true(size(fk(i).latlon,1),1);
    else
        lonidx=fk(i).latlon>=lonrng(1) & fk(i).latlon<=lonrng(2);
    end
    llidx=find(latidx & lonidx);
    
    % subsection
    sub=fk(i).beam(llidx,sidx,fidx);
    szsub=size(sub);
    
    % stats
    meddb.db(i)=median(sub(:))+fk(i).normdb;
    [mindb.db(i),idx]=min(sub(:));
    mindb.db(i)=mindb.db(i)+fk(i).normdb;
    [r,c,p]=ind2sub(szsub,idx);
    mindb.horzslow(i)=fk(i).horzslow(sidx(c));
    mindb.latlon(i,:)=fk(i).latlon(llidx(r),:);
    mindb.freq(i)=fk(i).freq(fidx(p));
    [maxdb.db(i),idx]=max(sub(:));
    maxdb.db(i)=maxdb.db(i)+fk(i).normdb;
    [r,c,p]=ind2sub(szsub,idx);
    maxdb.horzslow(i)=fk(i).horzslow(sidx(c));
    maxdb.latlon(i,:)=fk(i).latlon(llidx(r),:);
    maxdb.freq(i)=fk(i).freq(fidx(p));
end

end
