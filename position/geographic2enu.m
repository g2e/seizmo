function [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0,ellipsoid)
%GEOGRAPHIC2ENU    Converts geographic to local East/North/Up system
%
%    Usage:    [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0)
%              [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0,[a f])
%
%    Description: [E,N,U]=GEOGRAPHIC2ENU(LAT,LON,DEPTH,LAT0,LON0,DEPTH0)
%     converts coordinates in geographic latitude, longitude, depth to
%     local East/North/Up.  LAT0, LON0, DEPTH0 gives the reference point in
%     the local system (ie the origin point).  All latitudes & longitudes
%     must be in degrees.  Depths, E, N, U must be and are in kilometers.
%
%     [E,N,U]=GEOGRAPHIC2ENU(LAT,LON,DEPTH,LAT0,LON0,DEPTH0,[A F]) allows
%     specifying the ellipsoid parameters A (equatorial radius in
%     kilometers) and F (flattening).  The default corresponds to WGS-84.
%     This is compatible with output from Matlab's Mapping Toolbox function
%     ALMANAC.
%
%    Notes:
%
%    Examples:
%     Get local coordinates for a SEIZMO dataset:
%      st=getheader(data,'st');
%      [clat,clon]=arraycenter(st(:,1),st(:,2));
%      [e,n,u]=geographic2enu(st(:,1),st(:,2),(st(:,4)-st(:,3))/1000,...
%                             clat,clon,0);
%
%    See also: ENU2GEOGRAPHIC, GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC,
%              GEOCENTRIC2XYZ, XYZ2GEOCENTRIC,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        May   3, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   3, 2010 at 18:50 GMT

% todo:

% check nargin
msg=nargchk(6,7,nargin);
if(~isempty(msg)); error(msg); end

% default ellipsoid
if(nargin==6 || isempty(ellipsoid)); ellipsoid=[]; end

% check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(depth) ...
        || ~isreal(lat0) || ~isreal(lon0) || ~isreal(depth0))
    error('seizmo:geographic2enu:badInput',...
        'All position inputs must be real arrays!');
elseif(~isequalsizeorscalar(lat,lon,depth,lat0,lon0,depth0))
    error('seizmo:geographic2enu:badInput',...
        'All position inputs must be scalar or equal sized!');
end

% convert to xyz
[x,y,z]=geographic2xyz(lat,lon,depth,ellipsoid);
[x0,y0,z0]=geographic2xyz(lat0,lon0,depth0,ellipsoid);

% degrees to radians
d2r=pi/180;
lon0=lon0*d2r;
lat0=lat0*d2r;

% get local coordinates
e=-sin(lon0).*(x-x0)...
    +cos(lon0).*(y-y0);
n=-sin(lat0).*cos(lon0).*(x-x0)...
    -sin(lat0).*sin(lon0).*(y-y0)...
    +cos(lat0).*(z-z0);
u=cos(lat0).*cos(lon0).*(x-x0)...
    +cos(lat0).*sin(lon0).*(y-y0)...
    +sin(lat0).*(z-z0);

end
