function [lat]=geocentric2geographiclat(lat,ecc)
%GEOCENTRIC2GEOGRAPHICLAT    Convert latitude from geocentric to geographic
%
%    Usage:    lat=geocentric2geographiclat(lat)
%              lat=geocentric2geographiclat(lat,ecc)
%
%    Description:
%     LAT=GEOCENTRIC2GEOGRAPHICLAT(LAT) converts geocentric latitudes LAT
%     to geographic latitudes.  LAT is in degrees.  Assumes the WGS-84
%     reference ellipsoid.
%
%     LAT=GEOCENTRIC2GEOGRAPHICLAT(LAT,ECC) specifies the eccentricity for
%     the ellipsoid to use in the conversion.
%
%    Notes:
%     - If the location is not on the surface use GEOCENTRIC2GEOGRAPHIC.
%
%    Examples:
%     % Show the difference in latitudes (geographic pushes to the poles):
%     figure;
%     plot(-90:90,geocentric2geographiclat(-90:90)-(-90:90))
%     xlabel('geocentric latitude (^o)')
%     ylabel('geographic adjustment (^o)')
%
%    See also: GEOGRAPHIC2GEOCENTRICLAT, GEOGRAPHICLAT2RADIUS,
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
    error('seizmo:geocentric2geographiclat:badInput',...
        'LAT must be a real-valued array!');
elseif(~isreal(ecc) || (~isscalar(ecc) && ~isequal(size(ecc),size(lat)))...
        || any(ecc(:)>=1 | ecc(:)<0))
    error('seizmo:geocentric2geographiclat:badInput',...
        'ECC must be real-valued with 0<=ECC<1 !');
end

% convert to geographic
lat=atan2(sind(lat),(1-ecc^2).*cosd(lat)).*(180/pi);

end
