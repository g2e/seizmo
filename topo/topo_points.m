function [topo]=topo_points(lat,lon,toponame)
%TOPO_POINTS    Returns topography data at the points indicated
%
%    Usage:    topo=topo_points(lat,lon)
%              topo=topo_points(lat,lon,toponame)
%
%    Description:
%     TOPO=TOPO_POINTS(LAT,LON) grabs the topography (or z-value) for the
%     positions given by LAT & LON.  LAT & LON must be in degrees and must
%     be scalar or equal-sized arrays.  By default, this returns SRTM30+
%     topography in meters at the points given.
%
%     TOPO=TOPO_POINTS(LAT,LON,TOPONAME) allows changing the dataset to
%     pull the topography (z-values) from.  TOPONAME must be a string and
%     currently the valid strings are:
%      'srtm30plus' - SRTM30+ topography (30 arc-second resolution)
%                     (http://topex.ucsd.edu/WWW_html/srtm30_plus.html)
%      'etopo1_bed' - ETOPO1 bedrock topography (1 arc-minute resolution)
%      'etopo1_ice' - ETOPO1 ice topography (1 arc-minute resolution)
%      'etopo1_thk' - ETOPO1 ice thickness (1 arc-minute resolution)
%                     (http://www.ngdc.noaa.gov/mgg/global/global.html)
%
%    Notes:
%
%    Examples:
%     % Plot up an equatorial profile:
%     plot(-180:0.1:180,topo_points(0,-180:0.1:180))
%
%     % Compare SRTM30+ and ETOPO1:
%     plot(-180:0.1:180,topo_points(0,-180:0.1:180),'r',...
%          -180:0.1:180,topo_points(0,-180:0.1:180,'etopo1_bed'),'k:')
%
%    See also: TOPO_REGION, TOPO_COLORMAP

%     Version History:
%        Feb. 16, 2010 - initial version
%        May  17, 2010 - minor doc touch
%        May  19, 2010 - return nothing when given nothing
%        May  20, 2010 - added scrollbar (cause it is slow)
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 15:05 GMT

% todo:
% - etopo2, etopo5, crust2

% check nargin
error(nargchk(2,3,nargin));

% valid topo list (add new topo here)
validtopo={'srtm30plus' 'etopo1_bed' 'etopo1_ice' 'etopo1_thk'};

% check inputs
if(nargin==2 || isempty(toponame)); toponame='srtm30plus'; end
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:topo_points:badInput',...
        'LAT/LON must be real-valued!');
elseif(~isscalar(lat) && ~isscalar(lon) ...
        && ~isequal(size(lat),size(lon)))
    error('seizmo:topo_points:badInput',...
        'LAT/LON must be scalar or equal sized!');
elseif(~any(strcmpi(toponame,validtopo)))
    error('seizmo:topo_points:badInput',...
        ['TOPONAME must be one of the following:\n' ...
        sprintf('%s ',validtopo{:})]);
end

% expand scalars
if(isscalar(lat)); lat=lat(ones(size(lon))); end
if(isscalar(lon)); lon=lon(ones(size(lat))); end
npts=numel(lat);

% quick return if no points
if(~npts); topo=[]; return; end

% take care of wrap-around
[lat,lon]=fixlatlon(lat,lon);

% retrieve topo info
s=load(toponame,'version','registration',...
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

% verbosity
verbose=seizmoverbose;
if(verbose)
    cnt=0;
    disp('Getting Topography at Location(s)');
    print_time_left(cnt,npts);
end

% loop over tiles
topo=nan(size(lat)); cnt=0;
for i=1:numel(tilenames)
    % load tile
    z=load(toponame,tilenames{i});
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
    
    % detail message
    if(verbose); cnt=cnt+sum(pts); print_time_left(cnt,npts); end
end

end
