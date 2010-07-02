function [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0,ellipsoid)
%ENU2GEOGRAPHIC    Converts local East/North/Up system to geographic
%
%    Usage:    [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0)
%              [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0,[a f])
%
%    Description: [LAT,LON,DEPTH]=ENU2GEOGRAPHIC(E,N,U,LAT0,LON0,DEPTH0)
%     converts coordinates in local East/North/Up (with reference point
%     given by geographic position LAT0, LON0, DEPTH0) to their geographic
%     equivalent.  E, N, U, DEPTH0 must be in kilometers.  LAT0 & LON0 must
%     be in degrees.
%
%     [LAT,LON,DEPTH]=ENU2GEOGRAPHIC(E,N,U,LAT0,LON0,DEPTH0,[A F])
%     specifies the ellipsoid parameters A (equatorial radius in
%     kilometers) and F (flattening).  The default corresponds to WGS-84.
%     This is compatible with output from Matlab's Mapping Toolbox function
%     ALMANAC.
%
%    Notes:
%
%    Examples:
%     What altitude would you be at if you took a straight line North from
%     Saint Louis along the tangent plane for 2000km?:
%      [lat,lon,depth]=enu2geographic(0,2000,0,38.649,-90.305,0);
%      altitude=-depth
%
%    See also: GEOGRAPHIC2ENU, GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC,
%              GEOCENTRIC2XYZ, XYZ2GEOCENTRIC,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        May   3, 2010 - initial version
%        June 28, 2010 - fix nargchk
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 28, 2010 at 12:50 GMT

% todo:

% check nargin
error(nargchk(6,7,nargin));

% default ellipsoid
if(nargin==6 || isempty(ellipsoid)); ellipsoid=[]; end

% check inputs
if(~isreal(e) || ~isreal(n) || ~isreal(u) ...
        || ~isreal(lat0) || ~isreal(lon0) || ~isreal(depth0))
    error('seizmo:enu2geographic:badInput',...
        'All position inputs must be real arrays!');
elseif(~isequalsizeorscalar(e,n,u,lat0,lon0,depth0))
    error('seizmo:enu2geographic:badInput',...
        'All position inputs must be scalar or equal sized!');
end

% convert reference point to xyz
[x0,y0,z0]=geographic2xyz(lat0,lon0,depth0,ellipsoid);

% degrees to radians
d2r=pi/180;
lon0=lon0*d2r;
lat0=lat0*d2r;

% convert local ENU to xyz
x=-sin(lon0).*e-sin(lat0).*cos(lon0).*n+cos(lat0).*cos(lon0).*u+x0;
y=cos(lon0).*e-sin(lat0).*sin(lon0).*n+cos(lat0).*sin(lon0).*u+y0;
z=cos(lat0).*n+sin(lat0).*u+z0;

% xyz to geographic
[lat,lon,depth]=xyz2geographic(x,y,z,ellipsoid);

end
