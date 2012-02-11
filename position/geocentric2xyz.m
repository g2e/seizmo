function [x,y,z]=geocentric2xyz(lat,lon,radius,r)
%GEOCENTRIC2XYZ    Converts coordinates from geocentric to cartesian
%
%    Usage:    [x,y,z]=geocentric2xyz(lat,lon,radius)
%              [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    Description:
%     [X,Y,Z]=GEOCENTRIC2XYZ(LAT,LON,RADIUS) converts coordinates in
%     geocentric latitude, longitude, radius to Earth-centered, Earth-Fixed
%     (ECEF).  LAT and LON are in degrees.  X, Y, Z will match the units of
%     RADIUS.
%
%     [X,Y,Z]=GEOCENTRIC2XYZ(LAT,LON,DEPTH,R) allows specifying the radius
%     R of the sphere and DEPTH rather than RADIUS.  X, Y, Z will match the
%     units of DEPTH and R (must have the same units).
%
%    Notes:
%     - The ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    Examples:
%     % Find out how far a position is from the equatorial plane (z):
%     [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    See also: XYZ2GEOCENTRIC, XYZ2GEOGRAPHIC, GEOGRAPHIC2XYZ,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update, code cleaning, expand as necessary
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin));

% size up inputs
sz{1}=size(lat); sz{2}=size(lon); sz{3}=size(radius);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3});
if(nargin==4); sz{4}=size(r); n(4)=prod(sz{4}); end

% basic check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(radius))
    error('seizmo:geocentric2xyz:nonNumeric',...
        'All inputs must be reals!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:geocentric2xyz:badSize',...
        'All inputs must be equal sized or scalar!');
end

% check input (converts radius to depth)
if(nargin==4)
    if(~isnumeric(r))
        error('seizmo:geocentric2xyz:badR',...
        'R input must be real-valued!');
    end
    radius=r-radius; 
end

% convert to xyz
z=radius.*sind(lat);
x=radius.*cosd(lon).*cosd(lat);
y=radius.*sind(lon).*cosd(lat);

% expand z to lon if scalar
if(isscalar(z)); z=z(ones(sz{2})); end

end
