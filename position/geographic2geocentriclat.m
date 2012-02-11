function [lat]=geographic2geocentriclat(lat,ecc)
%GEOGRAPHIC2GEOCENTRICLAT    Convert latitude from geographic to geocentric
%
%    Usage:    lat=geographic2geocentriclat(lat)
%              lat=geographic2geocentriclat(lat,ecc)
%
%    Description:
%     LAT=GEOGRAPHIC2GEOCENTRICLAT(LAT) converts geographic latitudes LAT
%     to geocentric latitudes.  LAT is in degrees.  Assumes the WGS-84
%     reference ellipsoid.
%
%     LAT=GEOGRAPHIC2GEOCENTRICLAT(LAT,ECC) specifies the eccentricity
%     for the ellipsoid to use in the conversion.
%
%    Notes:
%     - If the location is not on the surface use GEOGRAPHIC2GEOCENTRIC.
%
%    Examples:
%     % Get the geocentric latitude for St. Louis, MO USA:
%     lat=geographic2geocentriclat(38.649)
%
%     % Show the difference in latitudes (authalic pushes to the equator):
%     figure;
%     plot(-90:90,geographic2geocentriclat(-90:90)-(-90:90))
%     xlabel('geographic latitude (^o)')
%     ylabel('geocentric adjustment (^o)')
%
%    See also: GEOCENTRIC2GEOGRAPHICLAT, GEOGRAPHICLAT2RADIUS,
%              AUTHALIC2GEOGRAPHICLAT, GEOGRAPHIC2AUTHALICLAT

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - minor doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update, allow nonscalar ecc
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% require 1 or 2 inputs
error(nargchk(1,2,nargin));

% assume WGS-84 ellipsoid if no eccentricity given
if(nargin==1); ecc=8.181919084262149e-02; end

% check inputs
if(~isreal(lat))
    error('seizmo:geographic2geocentriclat:badInput',...
        'LAT must be a real-valued array!');
elseif(~isreal(ecc) || (~isscalar(ecc) && ~isequal(size(ecc),size(lat)))...
        || any(ecc(:)>=1 | ecc(:)<0))
    error('seizmo:geographic2geocentriclat:badInput',...
        'ECC must be real-valued with 0<=ECC<1 !');
end

% convert to geocentric
lat=atan2((1-ecc^2).*sind(lat),cosd(lat)).*(180/pi);

end
