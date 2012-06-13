function [s]=geofsssub(s,frng,srng,latrng,lonrng)
%GEOFSSSUB    Extracts a subspectra of a geofss struct
%
%    Usage:    s=geofsssub(s,frng,srng,latrng,lonrng)
%
%    Description:
%     S=GEOFSSSUB(S,FRNG,SRNG,LATRNG,LONRNG) reduces the geofss spectra in
%     S to the ranges given by FRNG, SRNG, LATRNG, & LONRNG and returns S
%     with all pertinent info updated.  S is a struct (see GEOFSS for
%     details).  FRNG specifies the frequency range as [FREQLOW FREQHIGH]
%     in Hz.  SRNG is the horizontal slowness range as [SLOWLOW SLOWHIGH]
%     in sec/deg.  LATRNG & LONRNG provide the parallels and meridians that
%     define a window in the geographic domain and they are expected to be
%     in degrees as [LOW HIGH].  Note that longitudes positions will be
%     wrapped so all possible points are included in the LONRNG.
%
%     By default all ranges are set to [] (full range).
%
%    Notes:
%     - Although this will wrap the longitudes, it will NOT rearrange the
%       points to correspond to a monotonic grid required by PLOTGEOFSS.
%
%    Examples:
%     % Split a spectra into 10s period windows:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     s=geofss(d,[lat(:) lon(:)],27:33,[1/50 1/20],'center');
%     s1=geofsssub(s,[1/50 1/40]);
%     s2=geofsssub(s,[1/40 1/30]);
%     s3=geofsssub(s,[1/30 1/20]);
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSAVG, PLOTGEOFSS, GEOFSSFREQSLIDE,
%              GEOFSSSLOWSLIDE, GEOFSSFRAMESLIDE

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  8, 2012 - adapted from geofksubvol
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  8, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check struct
error(chkgeofss(s));

% don't allow averaged spectrums
if(any(~sum(reshape([s.vector],[2 numel(s)]))))
    error('seizmo:geofsssub:badInput',...
        'Cannot extract a subspectrum if S is already averaged!');
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
    error('seizmo:geofsssub:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofsssub:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isempty(latrng) && (~isreal(latrng) || numel(sla)~=2 || sla(1)~=1 ...
        || sla(2)~=2 || any(abs(latrng(:))>90)))
    error('seizmo:geofsssub:badInput',...
        'LATRNG must be a 1x2 array of [LATLOW LATHIGH] in degrees!');
end
if(~isempty(lonrng) && (~isreal(lonrng) || numel(slo)~=2 || slo(1)~=1 ...
        || slo(2)~=2))
    error('seizmo:geofsssub:badInput',...
        'LONRNG must be a 1x2 array of [LONLOW LONHIGH] in degrees!');
end

% extract subspectra
for i=1:numel(s)
    % freq range
    if(isempty(frng))
        fmin=min(s(i).freq);
        fmax=max(s(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    fidx=s(i).freq>=fmin & s(i).freq<=fmax;
    
    % slow range
    if(isempty(srng))
        smin=min(s(i).slow);
        smax=max(s(i).slow);
    else
        smin=srng(1);
        smax=srng(2);
    end
    sidx=s(i).slow>=smin & s(i).slow<=smax;
    
    % lat range
    if(isempty(latrng))
        latmin=min(s(i).latlon(:,1));
        latmax=max(s(i).latlon(:,1));
    else
        latmin=latrng(1);
        latmax=latrng(2);
    end
    latidx=s(i).latlon(:,1)>=latmin & s(i).latlon(:,1)<=latmax;
    
    % lon range
    if(isempty(lonrng))
        lonmin=min(s(i).latlon(:,2));
        lonmax=max(s(i).latlon(:,2));
    else
        lonmin=lonrng(1);
        lonmax=lonrng(2);
    end
    while(any(lonmin-s(i).latlon(:,2)>180))
        moveme=lonmin-s(i).latlon(:,2)>180;
        s(i).latlon(moveme,2)=s(i).latlon(moveme,2)+360;
    end
    while(any(s(i).latlon(:,2)-lonmax>180))
        moveme=s(i).latlon(:,2)-lonmax>180;
        s(i).latlon(moveme,2)=s(i).latlon(moveme,2)-360;
    end
    lonidx=s(i).latlon(:,2)>=lonmin & s(i).latlon(:,2)<=lonmax;
    
    % extract subspectra
    s(i).spectra=s(i).spectra(latidx & lonidx,sidx,fidx);
    
    % update info
    s(i).freq=s(i).freq(fidx);
    s(i).slow=s(i).slow(sidx);
    s(i).latlon=s(i).latlon(latidx & lonidx,:);
    if(isscalar(s(i).freq)); s(i).vector(1)=false; end
    if(isscalar(s(i).slow)); s(i).vector(2)=false; end
end

end
