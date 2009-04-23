function [lat,lon,depth]=geocentric2geodetic(lat,lon,depth,ellipsoid)
%GEOCENTRIC2GEODETIC    Converts coordinates from geocentric to geodetic
%
%    Usage:    [lat,lon,depth]=geocentric2geodetic(lat,lon,depth)
%              [lat,lon,depth]=geocentric2geodetic(lat,lon,depth,[a f])
%
%    Description: [LAT,LON,DEPTH]=GEOCENTRIC2GEODETIC(LAT,LON,DEPTH) 
%     converts arrays of coordinates from geocentric to geodetic latitude,
%     longitude, depth.  LAT and LON are in degrees.  DEPTH is in
%     kilometers.  The reference ellipsoid is assumed to be WGS-84.  The
%     volumetric radius is derived from ellipsoid parameters for the
%     spherical radius.
%
%     [LAT,LON,DEPTH]=GEOCENTRIC2GEODETIC(LAT,LON,DEPTH,[A F]) allows
%     specifying the ellipsoid parameters A (equatorial radius in 
%     kilometers) and F (flattening).  This is compatible with Matlab's 
%     Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     Get your earthquake location into geodetic coordinates:
%      [lat,lon,depth]=geocentric2geodetic(evla,evlo,evdp/1000)
%
%    See also: geodetic2geocentric, geocentric2xyz, xyz2geodetic

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:45 GMT

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
        error('seizmo:geocentric2geodetic:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get volumetric radius for this ellipsoid
r=(a^3*(1-f))^(1/3);

% convert to geodetic (via xyz)
[x,y,z]=geocentric2ecef(lat,lon,depth,r);
[lat,lon,depth]=ecef2geodetic(x,y,z,[a f]);

end
