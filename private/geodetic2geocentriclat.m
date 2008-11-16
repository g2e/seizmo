function [lat]=geodetic2geocentriclat(lat,ecc)
%GEODETIC2GEOCENTRICLAT    Convert latitudes from geodetic to geocentric
%
%    Description: GEODETIC2GEOCENTRICLAT(LATITUDES) converts LATITUDES that
%     are geodetic latitudes to geocentric latitudes.  LATITUDES units are
%     in degrees.  Assumes the WGS-84 reference ellipsoid.
%
%     GEODETIC2GEOCENTRICLAT(LATITUDES,ECC) specifies the eccentricity for
%     the ellipsoid to use in the conversion.
%
%    Notes:
%     - If the location is not on the surface use GEODETIC2GEOCENTRIC.
%
%    Tested on: Matlab r2007b
%
%    Usage:    latitudes=geodetic2geocentriclat(latitudes)
%              latitudes=geodetic2geocentriclat(latitudes,ecc)
%
%    Examples:
%     Get the geocentric latitude for St. Louis, MO USA:
%      latitude=geodetic2geocentriclat(38.649)
%
%    See also: geocentric2geodeticlat, geodeticlat2radius

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 10, 2008 at 08:40 GMT

% todo:

% require 1 or 2 inputs
error(nargchk(1,2,nargin))

% assume WGS-84 ellipsoid if no eccentricity given
if(nargin==1); ecc=8.181919084262149e-02; end

% check inputs
if(~isnumeric(lat) || isempty(lat))
    error('seizmo:geodetic2geocenticlat:badInput',...
        'LAT must be nonempty numeric array!');
elseif(~isnumeric(ecc) || ~isscalar(ecc) || ecc>=1 || ecc<0)
    error('seizmo:geodetic2geocentriclat:badInput',...
        'ECC must be numeric scalar with 0<=ECC<1 !');
end

% convert to geocentric
lat=atan2((1-ecc^2).*sind(lat),cosd(lat)).*(180/pi);

end
