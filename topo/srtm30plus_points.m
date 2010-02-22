function [topo]=srtm30plus_points(lat,lon)
%SRTM30PLUS_POINTS    Returns SRTM30plus data at the points indicated

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check inputs
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:points:badInput',...
        'LAT/LON must be real-valued!');
elseif(~isscalar(lat) && ~isscalar(lon) ...
        && ~isequal(size(lat),size(lon)))
    error('seizmo:points:badInput',...
        'LAT/LON must be scalar or equal sized!');
end

% expand scalars
if(isscalar(lat)); lat=lat(ones(size(lon))); end
if(isscalar(lon)); lon=lon(ones(size(lat))); end
npts=numel(lat);

% take care of wrap-around
[lat,lon]=fixlatlon(lat,lon);

% retrieve topo info
s=load('srtm30plus','version','registration',...
    'pixelsperdegree','latname','lonname');
ppd=s.pixelsperdegree;
spacing=1/ppd;

% get tile indices for each point
ll=90:-10:-90;
latidx=nan(size(lat));
for i=1:npts
    latidx(i)=find(lat(i)<=ll,1,'last');
end
latidx(latidx>18)=18;
ll=-180:10:180;
lonidx=nan(size(lat));
for i=1:npts
    lonidx(i)=find(lon(i)>=ll,1,'last');
end
lonidx(lonidx>36)=36;

% get tile names
[tidx,idx,idx]=unique([lonidx(:) latidx(:)],'rows');
tilenames=strcat(s.lonname(tidx(:,1))',s.latname(tidx(:,2))');

% loop over tiles
topo=nan(size(lat));
for i=1:numel(tilenames)
    % load tile
    z=load('srtm30plus',tilenames{i});
    z=z.(char(fieldnames(z)));
    
    % get lat/lon for this tile
    switch lower(s.registration)
        case 'pixel'
            % have to expand to assure all points are within tile
            z=[z(:,1) z z(:,end)]; %#ok
            z=[z(1,:); z; z(end,:)]; %#ok
            zlat=90-(tidx(i,2)-1)*10-spacing/2-(0:10*ppd-1)*spacing;
            zlon=-180+(tidx(i,1)-1)*10+spacing/2+(0:10*ppd-1)*spacing;
            zlat=[90-(tidx(i,2)-1)*10 zlat 90-(tidx(i,2))*10]; %#ok
            zlon=[-180+(tidx(i,1)-1)*10 zlon -180+(tidx(i,1))*10]; %#ok
        case 'grid'
            zlat=90-(tidx(i,2)-1)*10-(0:10*ppd)*spacing;
            zlon=-180+(tidx(i,1)-1)*10+(0:10*ppd)*spacing;
    end
    
    % points in this tile
    pts=idx==i;
    
    % get values
    topo(pts)=interp2(zlon,zlat,z,lon(pts),lat(pts),'nearest');
end

end
