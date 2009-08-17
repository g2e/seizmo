function [lat,lon,depth]=geodetic2geocentric(lat,lon,depth,ellipsoid)
%GEODETIC2GEOCENTRIC    Converts coordinates from geodetic to geocentric
%
%    Usage:    [lat,lon,depth]=geodetic2geocentric(lat,lon,depth)
%              [lat,lon,depth]=geodetic2geocentric(lat,lon,depth,[a f])
%
%    Description: [LAT,LON,DEPTH]=GEODETIC2GEOCENTRIC(LAT,LON,DEPTH) 
%     converts arrays of coordinates from geodetic to geocentric latitude,
%     longitude, depth.  LAT and LON are in degrees.  DEPTH is in
%     kilometers.  The reference ellipsoid is assumed to be WGS-84.  The
%     volumetric radius is derived from the ellipsoid parameters to find
%     the spherical depth.
%
%     [LAT,LON,DEPTH]=GEODETIC2GEOCENTRIC(LAT,LON,DEPTH,[A F]) allows
%     specifying the ellipsoid parameters A (equatorial radius in 
%     kilometers) and F (flattening).  This is compatible with Matlab's 
%     Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     Get your station locations into geocentric coordinates:
%      [lat,lon,depth]=geodetic2geocentric(stla,stlo,(stdp-stel)/1000)
%
%    See also: geocentric2geodetic, geodetic2xyz, xyz2geocentric

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% require 3 or 4 inputs
msg=nargchk(3,4,nargin);
if(~isempty(msg)); error(msg); end

% default - WGS-84 Reference Ellipsoid
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
        error('seizmo:geodetic2geocentric:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get volumetric radius for this ellipsoid
r=(a^3*(1-f))^(1/3);

% convert to geocentric (via xyz)
[x,y,z]=geodetic2xyz(lat,lon,depth,[a f]);
[lat,lon,depth]=xyz2geocentric(x,y,z,r);

end
