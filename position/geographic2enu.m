function [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0,ellipsoid)
%GEOGRAPHIC2ENU    Converts geographic to local East/North/Up system
%
%    Usage:    [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0)
%              [e,n,u]=geographic2enu(lat,lon,depth,lat0,lon0,depth0,[a f])
%
%    Description:
%     [E,N,U]=GEOGRAPHIC2ENU(LAT,LON,DEPTH,LAT0,LON0,DEPTH0) converts
%     coordinates in geographic latitude, longitude, depth to local
%     East/North/Up.  LAT0, LON0, DEPTH0 gives the reference point in the
%     local system (ie the origin point).  All latitudes & longitudes must
%     be in degrees.  E, N, U & depths are in kilometers.
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
%     % Get local coordinates for a SEIZMO dataset:
%     st=getheader(data,'st');
%     [clat,clon]=arraycenter(st(:,1),st(:,2));
%     [e,n,u]=geographic2enu(st(:,1),st(:,2),(st(:,4)-st(:,3))/1000,...
%                            clat,clon,0);
%
%    See also: ENU2GEOGRAPHIC, GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC,
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
sz{1}=size(lat); sz{2}=size(lon); sz{3}=size(depth);
sz{4}=size(lat0); sz{5}=size(lon0); sz{6}=size(depth0);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3});
n(4)=prod(sz{4}); n(5)=prod(sz{5}); n(6)=prod(sz{6});

% check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(depth) ...
        || ~isreal(lat0) || ~isreal(lon0) || ~isreal(depth0))
    error('seizmo:geographic2enu:badInput',...
        'All position inputs must be real arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:geographic2enu:badSize',...
        'All inputs must be equal sized or scalar!');
end

% convert to xyz
[x,y,z]=geographic2xyz(lat,lon,depth,ellipsoid);
[x0,y0,z0]=geographic2xyz(lat0,lon0,depth0,ellipsoid);

% for efficiency
slo=sind(lon0);
clo=cosd(lon0);
sla=sind(lat0);
cla=cosd(lat0);

% get local coordinates
e=-slo.*(x-x0)+clo.*(y-y0);
n=-sla.*clo.*(x-x0)-sla.*slo.*(y-y0)+cla.*(z-z0);
u=cla.*clo.*(x-x0)+cla.*slo.*(y-y0)+sla.*(z-z0);

end
