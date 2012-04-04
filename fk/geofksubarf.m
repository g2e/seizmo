function [arf]=geofksubarf(arf,srng,latrng,lonrng)
%GEOFKSUBARF    Extracts a subARF of a geofkarf volume
%
%    Usage:    arf=geofksubarf(arf,srng,latrng,lonrng)
%
%    Description:
%     ARF=GEOFKSUBARF(ARF,SRNG,LATRNG,LONRNG) extracts the geofkarf data in
%     ARF for the ranges given by FRNG, SRNG, LATRNG, & LONRNG and returns
%     the reduced array response function with all pertinent info updated.
%     ARF is a geofkarf struct (for struct details see that function).
%     SRNG defines the horizontal slowness range in seconds per degree as
%     [SLOWLOW SLOWHIGH].  LATRNG & LONRNG provide the parallels and
%     meridians that define a window in the geographic domain and they are
%     expected to be in degrees as [LOW HIGH].  Note that longitudes in the
%     response positions will be wrapped so all possible points are
%     included in the LONRNG.
%
%     By default all ranges are set to [] (full range).
%
%    Notes:
%     - Although this will wrap the longitudes, it will NOT rearrange the
%       points to correspond to a monotonic grid required by PLOTGEOFKARF.
%
%    Examples:
%
%    See also: GEOFKARF, PLOTGEOFKARF, GEOFKARF2MAP, GEOFKARFSLOWSLIDE,
%              UPDATEGEOFKARF, CHKGEOFKARFSTRUCT

%     Version History:
%        July  8, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:45 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow map
if(any(~sum(reshape([arf.volume],[2 numel(arf)]))))
    error('seizmo:geofksubarf:badInput',...
        'ARF must be a geofkarf volume!');
end

% check frequency range
if(nargin<2); srng=[]; end
if(nargin<3); latrng=[]; end
if(nargin<4); lonrng=[]; end
ss=size(srng);
sla=size(latrng);
slo=size(lonrng);
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofksubarf:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isempty(latrng) && (~isreal(latrng) || numel(sla)~=2 || sla(1)~=1 ...
        || sla(2)~=2 || any(abs(latrng(:))>90)))
    error('seizmo:geofksubarf:badInput',...
        'LATRNG must be a 1x2 array of [LATLOW LATHIGH] in degrees!');
end
if(~isempty(lonrng) && (~isreal(lonrng) || numel(slo)~=2 || slo(1)~=1 ...
        || slo(2)~=2))
    error('seizmo:geofksubarf:badInput',...
        'LONRNG must be a 1x2 array of [LONLOW LONHIGH] in degrees!');
end

% extract arf
for i=1:numel(arf)
    % slow range
    if(isempty(srng))
        smin=min(arf(i).horzslow);
        smax=max(arf(i).horzslow);
    else
        if(~arf(i).volume(1))
            error('seizmo:geofksubarf:badVol',...
                'Can not extract a slowness range from averaged volume!');
        end
        smin=srng(1);
        smax=srng(2);
    end
    if(arf(i).volume(1))
        sidx=arf(i).horzslow>=smin & arf(i).horzslow<=smax;
    else
        sidx=1;
    end
    
    % lat range
    if(isempty(latrng))
        latmin=min(arf(i).latlon(:,1));
        latmax=max(arf(i).latlon(:,1));
    else
        latmin=latrng(1);
        latmax=latrng(2);
    end
    latidx=arf(i).latlon(:,1)>=latmin & arf(i).latlon(:,1)<=latmax;
    
    % lon range
    if(isempty(lonrng))
        lonmin=min(arf(i).latlon(:,2));
        lonmax=max(arf(i).latlon(:,2));
    else
        lonmin=lonrng(1);
        lonmax=lonrng(2);
    end
    while(any(lonmin-arf(i).latlon(:,2)>180))
        moveme=lonmin-arf(i).latlon(:,2)>180;
        arf(i).latlon(moveme,2)=arf(i).latlon(moveme,2)+360;
    end
    while(any(arf(i).latlon(:,2)-lonmax>180))
        moveme=arf(i).latlon(:,2)-lonmax>180;
        arf(i).latlon(moveme,2)=arf(i).latlon(moveme,2)-360;
    end
    lonidx=arf(i).latlon(:,2)>=lonmin & arf(i).latlon(:,2)<=lonmax;
    
    % extract subarf
    arf(i).beam=arf(i).beam(latidx & lonidx,sidx);
    
    % update info
    maxdb=max(arf(i).beam(:));
    arf(i).beam=arf(i).beam-maxdb;
    arf(i).normdb=arf(i).normdb+maxdb;
    if(arf(i).volume(1))
        arf(i).horzslow=arf(i).horzslow(sidx);
    end
    arf(i).latlon=arf(i).latlon(latidx & lonidx,:);
end

end
