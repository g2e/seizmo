function [lat]=geocentric2geodeticlat(lat,ecc)
%GEOCENTRIC2GEODETICLAT    Convert latitudes from geocentric to geodetic
%
%    Usage:    latitudes=geocentric2geodeticlat(latitudes)
%              latitudes=geocentric2geodeticlat(latitudes,ecc)
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
%    Examples:
%     Show the difference in latitudes:
%      plot(geocentric2geodeticlat(-90:90)-(-90:90))
%
%    See also: geodetic2geocentriclat, geodeticlat2radius

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - minor doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% require 1 or 2 inputs
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% assume WGS-84 ellipsoid if no eccentricity given
if(nargin==1); ecc=8.181919084262149e-02; end

% check inputs
if(~isnumeric(lat) || isempty(lat))
    error('seizmo:geocentric2geodeticlat:badInput',...
        'LAT must be nonempty numeric array!');
elseif(~isnumeric(ecc) || ~isscalar(ecc) || ecc>=1 || ecc<0)
    error('seizmo:geocentric2geodeticlat:badInput',...
        'ECC must be numeric scalar with 0<=ECC<1 !');
end

% convert to geodetic
lat=atan2(sind(lat),(1-ecc^2).*cosd(lat)).*(180/pi);

end
