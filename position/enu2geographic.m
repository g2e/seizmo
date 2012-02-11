function [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0,ellipsoid)
%ENU2GEOGRAPHIC    Converts local East/North/Up system to geographic
%
%    Usage:    [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0)
%              [lat,lon,depth]=enu2geographic(e,n,u,lat0,lon0,depth0,[a f])
%
%    Description:
%     [LAT,LON,DEPTH]=ENU2GEOGRAPHIC(E,N,U,LAT0,LON0,DEPTH0) converts
%     coordinates in local East/North/Up (with reference point given by
%     geographic position LAT0, LON0, DEPTH0) to their geographic
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
%     % What altitude would you be at if you took a straight line North
%     % from Saint Louis along the tangent plane for 2000km?:
%     [lat,lon,depth]=enu2geographic(0,2000,0,38.649,-90.305,0);
%     altitude=-depth
%
%    See also: GEOGRAPHIC2ENU, GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC,
%              GEOCENTRIC2XYZ, XYZ2GEOCENTRIC,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        May   3, 2010 - initial version
%        June 28, 2010 - fix nargchk
%        Feb. 10, 2012 - doc update, code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 12:50 GMT

% todo:

% check nargin
error(nargchk(6,7,nargin));

% default ellipsoid
if(nargin==6 || isempty(ellipsoid)); ellipsoid=[]; end

% size up inputs
sz{1}=size(e); sz{2}=size(n); sz{3}=size(u);
sz{4}=size(lat0); sz{5}=size(lon0); sz{6}=size(depth0);
nm(1)=prod(sz{1}); nm(2)=prod(sz{2}); nm(3)=prod(sz{3});
nm(4)=prod(sz{4}); nm(5)=prod(sz{5}); nm(6)=prod(sz{6});

% check inputs
if(~isreal(e) || ~isreal(n) || ~isreal(u) ...
        || ~isreal(lat0) || ~isreal(lon0) || ~isreal(depth0))
    error('seizmo:enu2geographic:badInput',...
        'All position inputs must be real arrays!');
elseif(sum(nm~=1)>1 && ~isequal(sz{nm~=1}))
    error('seizmo:enu2geographic:badSize',...
        'All inputs must be equal sized or scalar!');
end

% convert reference point to xyz
[x0,y0,z0]=geographic2xyz(lat0,lon0,depth0,ellipsoid);

% for efficiency
slo=sind(lon0);
clo=cosd(lon0);
sla=sind(lat0);
cla=cosd(lat0);

% convert local ENU to xyz
x=-slo.*e-sla.*clo.*n+cla.*clo.*u+x0;
y=clo.*e-sla.*slo.*n+cla.*slo.*u+y0;
z=cla.*n+sla.*u+z0;

% xyz to geographic
[lat,lon,depth]=xyz2geographic(x,y,z,ellipsoid);

end
