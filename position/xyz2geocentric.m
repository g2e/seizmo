function [lat,lon,radius]=xyz2geocentric(x,y,z,r)
%XYZ2GEOCENTRIC    Converts coordinates from cartesian to geocentric
%
%    Usage:    [lat,lon,radius]=xyz2geocentric(x,y,z)
%              [lat,lon,depth]=xyz2geocentric(x,y,z,r)
%
%    Description:
%     [LAT,LON,RADIUS]=XYZ2GEOCENTRIC(X,Y,Z) converts arrays of coordinates
%     in Earth-centered, Earth-Fixed (ECEF) to geocentric latitude,
%     longitude, radius.  LAT and LON are in degrees.  X, Y and Z must have
%     the same units (so RADIUS will be in those units) and must be same
%     size arrays or scalars.
%
%     [LAT,LON,DEPTH]=XYZ2GEOCENTRIC(X,Y,Z,R) allows specifying the radius
%     R of the sphere so depth is returned rather than radius.  In this 
%     case, the units of X, Y, Z must match those of R (so DEPTH will be in
%     those units).
%
%    Notes:
%     - The ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    Examples:
%     % Find out how far an xyz position is from the earth's center:
%     [lat,lon,radius]=xyz2geocentric(x,y,z)
%
%    See also: GEOCENTRIC2XYZ, XYZ2GEOGRAPHIC, GEOGRAPHIC2XYZ,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Jan.  4, 2011 - minor code improvements
%        Feb. 10, 2012 - doc update, drop scalar r require, code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 20:15 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin));

% size up inputs
sz{1}=size(x); sz{2}=size(y); sz{3}=size(z);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3});
if(nargin==4); sz{4}=size(r); n(4)=prod(sz{4}); end

% basic check inputs
if(~isreal(x) || ~isreal(y) || ~isreal(z))
    error('seizmo:xyz2geocentric:nonReals','All inputs must be reals!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:xyz2geocentric:badSize',...
        'All inputs must be equal sized or scalar!');
end

% convert to geocentric
radius=sqrt(x.^2+y.^2+z.^2);
lon=atan2(y,x).*(180/pi);
lat=asind(z./radius);

% expand lon to z if scalar
if(isscalar(lon)); lon=lon(ones(sz{3})); end

% check input (converts radius to depth)
if(nargin==4)
    if(~isreal(r))
        error('seizmo:xyz2geocentric:badR',...
        'R input must be real-valued!');
    end
    radius=r-radius;
    
    % expand lat/lon to r if scalar
    if(isscalar(lat)); lat=lat(ones(sz{4})); end
    if(isscalar(lon)); lon=lon(ones(sz{4})); end
end

end
