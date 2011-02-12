function [lat]=authalic2geographiclat(lat,ecc)
%AUTHALIC2GEOGRAPHICLAT    Convert latitude from authalic to geographic
%
%    Usage:    latitudes=authalic2geographiclat(latitudes)
%              latitudes=authalic2geographiclat(latitudes,ecc)
%
%    Description: AUTHALIC2GEOGRAPHICLAT(LATITUDES) converts LATITUDES
%     that are authalic latitudes to geographic latitudes.  LATITUDES
%     units are in degrees.  Assumes the WGS-84 reference ellipsoid.
%
%     AUTHALIC2GEOGRAPHICLAT(LATITUDES,ECC) specifies the eccentricity
%     for the ellipsoid to use in the conversion.
%
%    Notes:
%
%    Examples:
%     Show the difference in latitudes:
%      plot(authalic2geographiclat(-90:90)-(-90:90))
%
%    See also: GEOGRAPHIC2AUTHALICLAT, GEOCENTRIC2GEOGRAPHICLAT,
%              GEOGRAPHIC2GEOCENTRICLAT

%     Version History:
%        Nov. 13, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% require 1 or 2 inputs
error(nargchk(1,2,nargin));

% assume WGS-84 ellipsoid if no eccentricity given
if(nargin==1); ecc=8.181919084262149e-02; end

% check inputs
if(~isnumeric(lat) || isempty(lat))
    error('seizmo:authalic2geographiclat:badInput',...
        'LAT must be nonempty numeric array!');
elseif(~isnumeric(ecc) || ~isscalar(ecc) || ecc>=1 || ecc<0)
    error('seizmo:authalic2geographiclat:badInput',...
        'ECC must be numeric scalar with 0<=ECC<1 !');
end

% convert to geographic
R2D=180/pi;
lat=lat/R2D;
lat=lat+((ecc^2)/3+31*(ecc^4)/180+517*(ecc^6)/5040)*sin(2*lat)...
    +(23*(ecc^4)/360+251*(ecc^6)/3780)*sin(4*lat)...
    +(761*(ecc^6)/45360)*sin(6*lat);
lat=lat*R2D;

end
