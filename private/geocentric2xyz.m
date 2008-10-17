function [x,y,z]=geocentric2xyz(lat,lon,radius,r)
%GEOCENTRIC2XYZ    Converts coordinates from geocentric to cartesian
%
%    Description: [X,Y,Z]=GEOCENTRIC2XYZ(LAT,LON,RADIUS) converts
%     coordinates in geocentric latitude, longitude, radius to 
%     Earth-centered, Earth-Fixed (ECEF).  LAT and LON are in degrees.  
%     X, Y, Z will match the units of RADIUS.
%
%     GEOCENTRIC2XYZ(LAT,LON,DEPTH,R) allows specifying the radius R of
%     the sphere and DEPTH rather than RADIUS.  X, Y, Z will match the
%     units of DEPTH and R (must have the same units).
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
%    Usage:    [x,y,z]=geocentric2xyz(lat,lon,radius)
%              [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    Examples:
%     Find out how far a position is from the equatorial plane (z):
%      [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    See also: xyz2geocentric, xyz2geodetic, geodetic2xyz

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 16:15 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin))

% check inputs
if(~isnumeric(lat) || ~isnumeric(lon) || ~isnumeric(radius))
    error('SAClab:geocentric2xyz:nonNumeric',...
    'All inputs must be numeric!');
elseif(isempty(lat) || ~isequal(size(lat),size(lon),size(radius)))
    error('SAClab:geocentric2xyz:unpairedCoord',...
        'Coordinate inputs must be nonempty, equal size arrays!');
end

% check input (converts depth to radius)
if(nargin==4)
    if(~isnumeric(r) || ~isscalar(r))
        error('SAClab:geocentric2xyz:badR',...
        'Sphere radius input must be numerical scalar!');
    end
    radius=r-radius; 
end

% convert to xyz
z=radius.*sind(lat);
x=radius.*cosd(lon).*cosd(lat);
y=radius.*sind(lon).*cosd(lat);

end
