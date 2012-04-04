function [topo,lat,lon]=topo_region(latrange,lonrange,toponame)
%TOPO_REGION    Returns topography data for the region indicated
%
%    Usage:    [topo,lat,lon]=topo_region(latrange,lonrange)
%              [topo,lat,lon]=topo_region(latrange,lonrange,toponame)
%
%    Description:
%     [TOPO,LAT,LON]=TOPO_REGION(LATRANGE,LONRANGE) returns a matrix of the
%     topography for the region bounded by LATRANGE & LONRANGE.  LATRANGE &
%     LONRANGE must in degrees and be 1x2 real-valued arrays.  If
%     LONRANGE(1)>LONRANGE(2) then the region wraps around the dateline.
%     The topography is from SRTM30+ and is in meters.
%
%     [TOPO,LAT,LON]=TOPO_REGION(LATRANGE,LONRANGE,TOPONAME) allows setting
%     the dataset that the topography (z-values) is pulled from.  TOPONAME
%     must be a string and the current possibilities are:
%      'srtm30plus' - SRTM30+ topography (30 arc-second resolution)
%                     (http://topex.ucsd.edu/WWW_html/srtm30_plus.html)
%      'etopo1_bed' - ETOPO1 bedrock topography (1 arc-minute resolution)
%      'etopo1_ice' - ETOPO1 ice topography (1 arc-minute resolution)
%      'etopo1_thk' - ETOPO1 ice thickness (1 arc-minute resolution)
%                     (http://www.ngdc.noaa.gov/mgg/global/global.html)
%    Notes:
%
%    Examples:
%     % The Bering Sea:
%     [topo,lat,lon]=topo_region([64 49],[162 -166],'etopo1_bed');
%     figure; surface(lon,lat,topo);
%     shading interp; topo_colormap(topo); axis equal tight;
%
%    See also: TOPO_POINTS, TOPO_COLORMAP

%     Version History:
%        Feb. 17, 2010 - initial version
%        May  17, 2010 - slight improvement to example
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% valid topo list (add new topo here)
validtopo={'srtm30plus' 'etopo1_bed' 'etopo1_ice' 'etopo1_thk'};

% check inputs
if(nargin==2 || isempty(toponame)); toponame='srtm30plus'; end
if(~isreal(latrange) || ~isreal(lonrange))
    error('seizmo:topo_region:badInput',...
        'LATRANGE/LONRANGE must be real-valued!');
elseif(~isequal(size(latrange),[1 2]) ...
        && ~isequal(size(lonrange),[1 2]))
    error('seizmo:topo_region:badInput',...
        'LATRANGE/LONRANGE must be 1x2!');
elseif(~any(strcmpi(toponame,validtopo)))
    error('seizmo:topo_region:badInput',...
        ['TOPONAME must be one of the following:\n' ...
        sprintf('%s ',validtopo{:})]);
end

% take care of wrap-around
[latrange,lonrange]=fixlatlon(latrange,lonrange);
latrange=sort(latrange);
lonflag=~isequal(lonrange,sort(lonrange));

% retrieve topo info
s=load(toponame,'version','registration','tilewidth',...
    'pixelsperdegree','latname','lonname');
tw=s.tilewidth;
ppd=s.pixelsperdegree;
spacing=1/ppd;

% what regions do we load?
% - always advance left to right
%   so we can go around dateline
lat1=90:-tw:-90+tw; lat2=90-tw:-tw:-90;
lon1=-180:tw:180-tw; lon2=-180+tw:tw:180;
inlat=find(lat1>latrange(1) & lat2<latrange(2));
if(lonflag)
    inlon=[find(lonrange(1)<lon2) find(lonrange(2)>lon1)];
else
    inlon=find(lon2>lonrange(1) & lon1<lonrange(2));
end

% loop over tiles
nlon=numel(inlon);
nlat=numel(inlat);
topo=cell(nlat,nlon);
for i=1:nlon
    for j=1:nlat
        % load tile
        topo{j,i}=load(toponame,[s.lonname{inlon(i)} s.latname{inlat(j)}]);
        topo{j,i}=topo{j,i}.(char(fieldnames(topo{j,i})));
        
        % trim edges of grid registered if not far left or bottom
        switch lower(s.registration)
            case 'grid'
                if(i~=nlon)
                    topo{j,i}=topo{j,i}(:,1:end-1);
                end
                if(j~=nlat)
                    topo{j,i}=topo{j,i}(1:end-1,:);
                end
        end
    end
end

% combine
topo=cell2mat(topo);

% get lat/lon
switch lower(s.registration)
    case 'pixel'
        lat=90-(inlat(1)-1)*tw-spacing/2-(0:nlat*tw*ppd-1)*spacing;
        lon=-180+(inlon(1)-1)*tw+spacing/2+(0:nlon*tw*ppd-1)*spacing;
    case 'grid'
        lat=90-(inlat(1)-1)*tw-(0:nlat*tw*ppd)*spacing;
        lon=-180+(inlon(1)-1)*tw+(0:nlon*tw*ppd)*spacing;
end

% truncate grid
[nlat,nlon]=size(topo);
inlat=[find(lat<latrange(2),1,'first')-1 ...
    find(lat>latrange(1),1,'last')+1];
inlat(inlat==0)=1; inlat(inlat>nlat)=nlat;
if(lonflag)
    inlon=[find(lon>lonrange(1),1,'first')-1 ...
        find(lon<(lonrange(2)+360),1,'last')+1];
    inlon(inlon==0)=1; inlon(inlon>nlon)=nlon;
else
    inlon=[find(lon>lonrange(1),1,'first')-1 ...
        find(lon<lonrange(2),1,'last')+1];
    inlon(inlon==0)=1; inlon(inlon>nlon)=nlon;
end
lat=lat(inlat(1):inlat(2));
lon=lon(inlon(1):inlon(2));
topo=double(topo(inlat(1):inlat(2),inlon(1):inlon(2)));

end
