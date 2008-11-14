function [lat]=geocentric2geodeticlat(lat,ecc)
%GEOCENTRIC2GEODETICLAT    Convert latitudes from geocentric to geodetic
%
%    Description: GEOCENTRIC2GEODETICLAT(LATITUDES) converts LATITUDES that
%     are geocentric latitudes to geodetic latitudes.  LATITUDES units are
%     in degrees.  Assumes the WGS-84 reference ellipsoid.
%
%     GEOCENTRIC2GEODETICLAT(LATITUDES,ECC) specifies the eccentricity for
%     the ellipsoid to use in the conversion.
%
%    Notes:
%     - If the location is not on the surface use GEOCENTRIC2GEODETIC.
%
%    Tested on: Matlab r2007b
%
%    Usage:    latitudes=geocentric2geodeticlat(latitudes)
%              latitudes=geocentric2geodeticlat(latitudes,ecc)
%
%    Examples:
%     Show the difference in latitudes:
%      plot(geocentric2geodeticlat(-90:90)-(-90:90))
%
%    See also: geodetic2geocentriclat, geodeticlat2radius

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
    error('SAClab:geocentric2geodeticlat:badInput',...
        'LAT must be nonempty numeric array!');
elseif(~isnumeric(ecc) || ~isscalar(ecc) || ecc>=1 || ecc<0)
    error('SAClab:geocentric2geodeticlat:badInput',...
        'ECC must be numeric scalar with 0<=ECC<1 !');
end

% convert to geodetic
lat=atan2(sind(lat),(1-ecc^2).*cosd(lat)).*(180/pi);

end
