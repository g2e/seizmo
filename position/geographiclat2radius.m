function [radius]=geographiclat2radius(lat,ellipsoid)
%GEOGRAPHICLAT2RADIUS    Returns the radius at a geographic latitude
%
%    Usage:    radius=geographiclat2radius(latitudes)
%              radius=geographiclat2radius(latitudes,[a f])
%
%    Description:
%     GEOGRAPHICLAT2RADIUS(LATITUDES) returns the radii at geographic
%     latitudes LATITUDES.  LATITUDES must be in degrees.  Assumes the
%     WGS-84 reference ellipsoid.
%
%     GEOGRAPHICLAT2RADIUS(LATITUDES,[A F]) allows specifying the ellipsoid
%     parameters A (equatorial radius in kilometers) and F (flattening).
%     This is compatible with Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     % Get the radius for St. Louis, MO USA:
%     radius=geographiclat2radius(38.649)
%
%    See also: GEOCENTRIC2GEOGRAPHICLAT, GEOGRAPHIC2GEOCENTRICLAT

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - minor doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% require 1 or 2 inputs
error(nargchk(1,2,nargin));

% default ellipsoid - WGS-84 Reference Ellipsoid
if(nargin==1 || isempty(ellipsoid))
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
        error('seizmo:geographiclat2radius:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get ellipsoid parameters
b=a*(1-f);
a2=a^2;
b2=b^2;

% optimization
b2sin2lat=b2.*sind(lat).^2;
a2cos2lat=a2.*cosd(lat).^2;

% get radius
radius=sqrt((a2.*a2cos2lat+b2.*b2sin2lat)./(a2cos2lat+b2sin2lat));

end
