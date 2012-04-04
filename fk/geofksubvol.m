function [vol]=geofksubvol(vol,frng,srng,latrng,lonrng)
%GEOFKSUBVOL    Extracts a subvolume of a geofk volume
%
%    Usage:    vol=geofksubvol(vol,frng,srng,latrng,lonrng)
%
%    Description:
%     VOL=GEOFKSUBVOL(VOL,FRNG,SRNG,LATRNG,LONRNG) extracts the geofk beam
%     data in VOL for the ranges given by FRNG, SRNG, LATRNG, & LONRNG and
%     returns the reduced volume with all pertinent info updated.  VOL is a
%     geofk struct (see GEOFKXCVOLUME for details).  FRNG gives the
%     frequency range to extract as [FREQLOW FREQHIGH] in Hz.  SRNG defines
%     the horizontal slowness range in seconds per degree as
%     [SLOWLOW SLOWHIGH].  LATRNG & LONRNG provide the parallels and
%     meridians that define a window in the geographic domain and they are
%     expected to be in degrees as [LOW HIGH].  Note that longitudes in the
%     beam data positions will be wrapped so all possible points are
%     included in the LONRNG.
%
%     By default all ranges are set to [] (full range).
%
%    Notes:
%     - Although this will wrap the longitudes, it will NOT rearrange the
%       points to correspond to a monotonic grid required by PLOTGEOFKMAP.
%
%    Examples:
%     % Split a volume into 10s period windows:
%     svol=fkvolume(data,50,201,[1/50 1/20]);
%     svol1=geofksubvol(svol,[1/50 1/40]);
%     svol2=geofksubvol(svol,[1/40 1/30]);
%     svol3=geofksubvol(svol,[1/30 1/20]);
%     svol4=geofksubvol(svol,[1/20 1/10]);
%
%    See also: GEOFKFREQSLIDE, GEOFKSLOWSLIDE, PLOTGEOFKMAP, GEOFKXCVOLUME,
%              GEOFKVOL2MAP, CHKGEOFKSTRUCT

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check fk struct
error(chkgeofkstruct(vol));

% don't allow map
if(any(~sum(reshape([vol.volume],[2 numel(vol)]))))
    error('seizmo:geofksubvol:badInput',...
        'VOL must be a geofk beam volume!');
end

% check frequency range
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
    error('seizmo:geofksubvol:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofksubvol:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isempty(latrng) && (~isreal(latrng) || numel(sla)~=2 || sla(1)~=1 ...
        || sla(2)~=2 || any(abs(latrng(:))>90)))
    error('seizmo:geofksubvol:badInput',...
        'LATRNG must be a 1x2 array of [LATLOW LATHIGH] in degrees!');
end
if(~isempty(lonrng) && (~isreal(lonrng) || numel(slo)~=2 || slo(1)~=1 ...
        || slo(2)~=2))
    error('seizmo:geofksubvol:badInput',...
        'LONRNG must be a 1x2 array of [LONLOW LONHIGH] in degrees!');
end

% extract volume
for i=1:numel(vol)
    % freq range
    if(isempty(frng))
        fmin=min(vol(i).freq);
        fmax=max(vol(i).freq);
    else
        if(~vol(i).volume(2))
            error('seizmo:geofksubvol:badVol',...
                'Can not extract a freq range from averaged volume!');
        end
        fmin=frng(1);
        fmax=frng(2);
    end
    if(vol(i).volume(2))
        fidx=vol(i).freq>=fmin & vol(i).freq<=fmax;
    else
        fidx=1;
    end
    
    % slow range
    if(isempty(srng))
        smin=min(vol(i).horzslow);
        smax=max(vol(i).horzslow);
    else
        if(~vol(i).volume(1))
            error('seizmo:geofksubvol:badVol',...
                'Can not extract a slowness range from averaged volume!');
        end
        smin=srng(1);
        smax=srng(2);
    end
    if(vol(i).volume(1))
        sidx=vol(i).horzslow>=smin & vol(i).horzslow<=smax;
    else
        sidx=1;
    end
    
    % lat range
    if(isempty(latrng))
        latmin=min(vol(i).latlon(:,1));
        latmax=max(vol(i).latlon(:,1));
    else
        latmin=latrng(1);
        latmax=latrng(2);
    end
    latidx=vol(i).latlon(:,1)>=latmin & vol(i).latlon(:,1)<=latmax;
    
    % lon range
    if(isempty(lonrng))
        lonmin=min(vol(i).latlon(:,2));
        lonmax=max(vol(i).latlon(:,2));
    else
        lonmin=lonrng(1);
        lonmax=lonrng(2);
    end
    while(any(lonmin-vol(i).latlon(:,2)>180))
        moveme=lonmin-vol(i).latlon(:,2)>180;
        vol(i).latlon(moveme,2)=vol(i).latlon(moveme,2)+360;
    end
    while(any(vol(i).latlon(:,2)-lonmax>180))
        moveme=vol(i).latlon(:,2)-lonmax>180;
        vol(i).latlon(moveme,2)=vol(i).latlon(moveme,2)-360;
    end
    lonidx=vol(i).latlon(:,2)>=lonmin & vol(i).latlon(:,2)<=lonmax;
    
    % extract subvolume
    vol(i).beam=vol(i).beam(latidx & lonidx,sidx,fidx);
    
    % update info
    maxdb=max(vol(i).beam(:));
    vol(i).beam=vol(i).beam-maxdb;
    vol(i).normdb=vol(i).normdb+maxdb;
    if(vol(i).volume(2))
        vol(i).freq=vol(i).freq(fidx);
    end
    if(vol(i).volume(1))
        vol(i).horzslow=vol(i).horzslow(sidx);
    end
    vol(i).latlon=vol(i).latlon(latidx & lonidx,:);
end

end
