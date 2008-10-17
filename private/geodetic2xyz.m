function [x,y,z]=geodetic2xyz(lat,lon,depth,ellipsoid)
%GEODETIC2XYZ    Converts coordinates from geodetic to cartesian
%
%    Description: [X,Y,Z]=GEODETIC2XYZ(LAT,LON,DEPTH) converts coordinates
%     in geodetic latitude, longitude, depth to Earth-centered, Earth-fixed
%     (ECEF).  LAT and LON are in degrees.  DEPTH, X, Y, Z must be/are in
%     kilometers.  The reference ellipsoid is assumed to be WGS-84.
%
%     GEODETIC2ECEF(LAT,LON,DEPTH,[A F]) allows specifying the ellipsoid
%     parameters A (equatorial radius in kilometers) and F (flattening).
%     This is compatible with Matlab's Mapping Toolbox function ALMANAC.
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
%    Usage:    [x,y,z]=geodetic2xyz(lat,lon,depth)
%              [x,y,z]=geodetic2xyz(lat,lon,depth,[a f])
%
%    Examples:
%     Get out of the geodetic coordinate system and into cartesian:
%      [x,y,z]=geodetic2xyz(lat,lon,depth)
%
%    See also: xyz2geodetic, spherical2xyz, xyz2spherical

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 16:35 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin))

% default - WGS-84 Reference Ellipsoid
if(nargin==3)
    % a=radius at equator (major axis)
    % f=flattening
    a=6378.137;
    f=1/298.257223563;
else
    % manually specify ellipsoid (will accept almanac output)
    if(isnumeric(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1);
        f=ellipsoid(2);
    else
        error('SAClab:geodetic2xyz:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% check inputs
if(~isnumeric(lat) || ~isnumeric(lon) || ~isnumeric(depth))
    error('SAClab:geodetic2xyz:nonNumeric','All inputs must be numeric!');
elseif(isempty(lat) || ~isequal(size(lat),size(lon),size(depth)))
    error('SAClab:geodetic2xyz:unpairedCoord',...
        'Coordinate inputs must be nonempty, equal size arrays!');
end

% optimization (reduce CPU USAGE vs MEMORY)
e2=f*(2-f);
sinlat=sind(lat);
achi=a./sqrt(1-e2.*sinlat.^2);
c=(achi-depth).*cosd(lat);

% convert to xyz
x=c.*cosd(lon);
y=c.*sind(lon);
z=(achi.*(1-e2)-depth).*sinlat;

end
