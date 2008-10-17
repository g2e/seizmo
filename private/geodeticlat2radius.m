function [radius]=geodeticlat2radius(lat,ellipsoid)
%GEODETICLAT2RADIUS    Returns the radius at a geodetic latitude
%
%    Description: GEODETICLAT2RADIUS(LATITUDES) returns the radii at
%     geodetic latitudes LATITUDES.  LATITUDES must be in degrees.  Assumes
%     the WGS-84 reference ellipsoid.
%
%     GEODETICLAT2RADIUS(LATITUDES,[A F]) allows specifying the ellipsoid
%     parameters A (equatorial radius in kilometers) and F (flattening).
%     This is compatible with Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    radius=geodeticlat2radius(latitudes)
%              radius=geodeticlat2radius(latitudes,[a f])
%
%    Examples:
%     Get the radius for St. Louis, MO USA:
%      radius=geodeticlat2radius(38.649)
%
%    See also: geocentric2geodeticlat, geodetic2geocentriclat

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 17:30 GMT

% todo:

% require 1 or 2 inputs
error(nargchk(1,2,nargin))

% default - WGS-84 Reference Ellipsoid
if(nargin==1)
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
        error('SAClab:geodeticlat2radius:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get ellipsoid parameters
b=a*(1-f);
a2=a^2;
b2=b^2;

% convert to radians
lat=lat.*(pi/180);

% optimization
b2sin2lat=b2.*sin(lat).^2;
a2cos2lat=a2.*cos(lat).^2;

% get radius
radius=sqrt((a2.*a2cos2lat+b2.*b2sin2lat)./(a2cos2lat+b2sin2lat));

end
