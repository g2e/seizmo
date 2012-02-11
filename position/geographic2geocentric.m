function [lat,lon,depth]=geographic2geocentric(lat,lon,depth,ellipsoid)
%GEOGRAPHIC2GEOCENTRIC    Converts coordinates from geographic 2 geocentric
%
%    Usage:    [lat,lon,depth]=geographic2geocentric(lat,lon,depth)
%              [lat,lon,depth]=geographic2geocentric(lat,lon,depth,[a f])
%
%    Description:
%     [LAT,LON,DEPTH]=GEOGRAPHIC2GEOCENTRIC(LAT,LON,DEPTH) converts arrays
%     of coordinates from geographic to geocentric latitude, longitude,
%     and depth.  LAT and LON are in degrees.  DEPTH is in kilometers.  The
%     reference ellipsoid is assumed to be WGS-84.  The volumetric radius
%     is derived from the ellipsoid parameters to find the spherical depth.
%
%     [LAT,LON,DEPTH]=GEOGRAPHIC2GEOCENTRIC(LAT,LON,DEPTH,[A F]) allows
%     specifying the ellipsoid parameters A (equatorial radius in 
%     kilometers) and F (flattening).  This is compatible with Matlab's 
%     Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     % Get your station locations into geocentric coordinates:
%     [lat,lon,depth]=geographic2geocentric(stla,stlo,(stdp-stel)/1000)
%
%    See also: GEOCENTRIC2GEOGRAPHIC, GEOGRAPHIC2XYZ, XYZ2GEOCENTRIC,
%              XYZ2GEOGRAPHIC, GEOCENTRIC2XYZ

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update
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
    if(isnumeric(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1);
        f=ellipsoid(2);
    else
        error('seizmo:geographic2geocentric:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get volumetric radius for this ellipsoid
r=(a^3*(1-f))^(1/3);

% convert to geocentric (via xyz)
[x,y,z]=geographic2xyz(lat,lon,depth,[a f]);
[lat,lon,depth]=xyz2geocentric(x,y,z,r);

end
