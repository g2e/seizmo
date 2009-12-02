function [lat]=geographic2authaliclat(lat,ecc)
%GEOGRAPHIC2AUTHALICLAT    Convert latitude from geographic to authalic
%
%    Usage:    latitudes=geographic2authaliclat(latitudes)
%              latitudes=geographic2authaliclat(latitudes,ecc)
%
%    Description: GEOGRAPHIC2AUTHALICLAT(LATITUDES) converts LATITUDES
%     that are geographic latitudes to authalic latitudes.  LATITUDES
%     units are in degrees.  Assumes the WGS-84 reference ellipsoid.
%
%     GEOGRAPHIC2AUTHALICLAT(LATITUDES,ECC) specifies the eccentricity
%     for the ellipsoid to use in the conversion.
%
%    Notes:
%
%    Examples:
%     Get the authalic latitude for St. Louis, MO USA:
%      latitude=geographic2authaliclat(38.649)
%
%    See also: AUTHALIC2GEOGRAPHICLAT, GEOGRAPHIC2GEOCENTRICLAT,
%              GEOCENTRIC2GEOGRAPHICLAT

%     Version History:
%        Nov. 13, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2009 at 20:55 GMT

% todo:

% require 1 or 2 inputs
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% assume WGS-84 ellipsoid if no eccentricity given
if(nargin==1); ecc=8.181919084262149e-02; end

% check inputs
if(~isnumeric(lat) || isempty(lat))
    error('seizmo:geographic2authaliclat:badInput',...
        'LAT must be nonempty numeric array!');
elseif(~isnumeric(ecc) || ~isscalar(ecc) || ecc>=1 || ecc<0)
    error('seizmo:geographic2authaliclat:badInput',...
        'ECC must be numeric scalar with 0<=ECC<1 !');
end

% convert to authalic
R2D=180/pi;
lat=lat/R2D;
lat=lat+((ecc^2)/3+31*(ecc^4)/180+59*(ecc^6)/560)*sin(2*lat)...
    +(17*(ecc^4)/360+61*(ecc^6)/1260)*sin(4*lat)...
    +(383*(ecc^6)/45360)*sin(6*lat);
lat=lat*R2D;

end
