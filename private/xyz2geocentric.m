function [lat,lon,radius]=xyz2geocentric(x,y,z,r)
%XYZ2SPHERICAL    Converts coordinates from cartesian to geocentric
%
%    Description: [LAT,LON,RADIUS]=XYZ2SPHERICAL(X,Y,Z) converts arrays of
%     coordinates in Earth-centered, Earth-Fixed (ECEF) to geocentric
%     latitude, longitude, radius.  LAT and LON are in degrees.  X, Y and Z
%     must be the same size and in the same units so that RADIUS will be
%     as well.
%
%     [LAT,LON,DEPTH]=XYZ2SPHERICAL(X,Y,Z,R) allows specifying the radius
%     R of the sphere so depth is returned rather than radius.  In this 
%     case, the size and units of X, Y, Z must match those of R so that 
%     DEPTH will match in size and units as well.
%
%    Notes:
%     - the ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    [lat,lon,radius]=xyz2geocentric(x,y,z)
%              [lat,lon,depth]=xyz2geocentric(x,y,z,r)
%
%    Examples:
%     Find out how far an xyz position is from the earth's center:
%      [lat,lon,radius]=xyz2geocentric(x,y,z)
%
%    See also: geocentric2xyz, xyz2geodetic, geodetic2xyz

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 16:45 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin))

% check inputs
if(~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z))
    error('SAClab:xyz2geocentric:nonNumeric',...
    'All inputs must be numeric!');
elseif(isempty(x) || ~isequal(size(x),size(y),size(z)))
    error('SAClab:xyz2geocentric:unpairedCoord',...
        'Coordinate inputs must be nonempty, equal size arrays!');
end

% convert to geocentric
radius=sqrt(x.^2+y.^2+z.^2);
lon=atan2(y,x).*(180/pi);
lat=asind(z./radius);

% check input (converts radius to depth)
if(nargin==4)
    if(~isnumeric(r) || ~isscalar(r))
        error('SAClab:xyz2geocentric:badR',...
        'R input must be numerical scalar!');
    end
    radius=r-radius; 
end

end
