function [lat,lon,depth]=xyz2geographic(x,y,z,ellipsoid)
%XYZ2GEOGRAPHIC    Converts coordinates from cartesian to geographic
%
%    Usage:    [lat,lon,depth]=xyz2geographic(x,y,z)
%              [lat,lon,depth]=xyz2geographic(x,y,z,[a f])
%
%    Description:
%     [LAT,LON,DEPTH]=XYZ2GEOGRAPHIC(X,Y,Z) converts coordinates in
%     Earth-centered, Earth-fixed (ECEF) coordinates to geographic
%     latitude, longitude, depth.  LAT and LON are in degrees.  DEPTH, X,
%     Y, Z must be/are in kilometers.  The reference ellipsoid is assumed
%     to be WGS-84.
%
%     [LAT,LON,DEPTH]=XYZ2GEOGRAPHIC(X,Y,Z,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  This is compatible with output from Matlab's Mapping
%     Toolbox function ALMANAC.  By default the ellipsoid parameters are
%     set to those of the reference ellipsoid WGS-84.
%
%    Notes:
%     - Utilizes the preferred algorithm in:
%        J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates to 
%        geographic coordinates," Aerospace and Electronic Systems, IEEE 
%        Transactions on, vol. 30, pp. 957-961, 1994
%     - The ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    Examples:
%     % Find the geographic position of some point given in xyz:
%     [lat,lon,depth]=xyz2geographic(3000,3000,3000)
%
%    See also: GEOGRAPHIC2XYZ, GEOCENTRIC2XYZ, XYZ2GEOCENTRIC,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - name change: geodetic to geographic
%        May   3, 2010 - better checking
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update, minor code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin));

% default ellipsoid - WGS-84 Reference Ellipsoid
if(nargin==3 || isempty(ellipsoid))
    % a=radius at equator (major axis)
    % f=flattening
    a=6378.137;
    f=1/298.257223563;
else
    % manually specify ellipsoid (will accept almanac output)
    if(isreal(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1);
        f=ellipsoid(2);
    else
        error('seizmo:xyz2geographic:badEllipsoid',...
            ['Ellipsoid must be a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% size up inputs
sz{1}=size(x); sz{2}=size(y); sz{3}=size(z);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3});

% basic check inputs
if(~isreal(x) || ~isreal(y) || ~isreal(z))
    error('seizmo:xyz2geographic:nonNumeric',...
        'All inputs must be numeric!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:xyz2geographic:badSize',...
        'All inputs must be equal sized or scalar!');
end

% vectorized setup
f12=(1-f)^2;
a2=a^2;
b2=a2*f12;
e2=f*(2-f);
z2=z.^2;
z21e2=z2.*f12;
r=sqrt(x.^2+y.^2);
r2=r.^2;
f=54.*b2.*z2;
g=r2+z21e2-e2.*(a2-b2);
c=(e2^2.*f.*r2)./(g.^3);
s=(1+c+sqrt(c.^2+2.*c)).^(1/3);
p=f./(3.*(s+1./s+1).^2.*g.^2); clear c s g f;
q=sqrt(1+2.*e2^2.*p);
q1=1./q;
q11=1./(1+q);
ro=-p.*e2.*r.*q11+sqrt(0.5.*a2.*(1+q1)-p.*z21e2.*q1.*q11-0.5.*p.*r2);
re2ro2=(r-e2.*ro).^2; clear ro r2 p q q1 q11;
u=sqrt(re2ro2+z2);
b2av=b2./(a.*sqrt(re2ro2+z21e2)); clear z2 re2ro2 z21e2;
r2d=180/pi;

% get geographic coords
depth=-u.*(1-b2av);
lat=atan((z+e2.*b2av.*z./f12)./r).*r2d;
lon=atan2(y,x).*r2d;

% expand lon to size of z if scalar
if(isscalar(lon)); lon=lon(ones(sz{3})); end

end
