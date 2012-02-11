function [lat]=authalic2geographiclat(lat,ecc)
%AUTHALIC2GEOGRAPHICLAT    Convert latitude from authalic to geographic
%
%    Usage:    latitudes=authalic2geographiclat(latitudes)
%              latitudes=authalic2geographiclat(latitudes,ecc)
%
%    Description:
%     AUTHALIC2GEOGRAPHICLAT(LATITUDES) converts LATITUDES that are
%     authalic latitudes to geographic latitudes.  LATITUDES units are in
%     degrees.  Assumes the WGS-84 reference ellipsoid.
%
%     AUTHALIC2GEOGRAPHICLAT(LATITUDES,ECC) specifies the eccentricity
%     for the ellipsoid to use in the conversion.
%
%    Notes:
%
%    Examples:
%     % Show the difference in latitudes (geographic pushes to the poles):
%     figure;
%     plot(-90:90,authalic2geographiclat(-90:90)-(-90:90))
%     xlabel('authalic latitude (^o)')
%     ylabel('geographic adjustment (^o)')
%
%    See also: GEOGRAPHIC2AUTHALICLAT, GEOCENTRIC2GEOGRAPHICLAT,
%              GEOGRAPHIC2GEOCENTRICLAT

%     Version History:
%        Nov. 13, 2009 - initial version
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
    error('seizmo:authalic2geographiclat:badInput',...
        'LAT must be a real-valued array!');
elseif(~isreal(ecc) || (~isscalar(ecc) && ~isequal(size(ecc),size(lat)))...
        || any(ecc(:)>=1 | ecc(:)<0))
    error('seizmo:authalic2geographiclat:badInput',...
        'ECC must be real-valued with 0<=ECC<1 !');
end

% convert to geographic
R2D=180/pi;
lat=lat+R2D*(((ecc^2)/3+31*(ecc^4)/180+517*(ecc^6)/5040)*sind(2*lat)...
    +(23*(ecc^4)/360+251*(ecc^6)/3780)*sind(4*lat)...
    +(761*(ecc^6)/45360)*sind(6*lat));

end
