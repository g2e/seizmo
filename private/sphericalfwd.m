function [stla,stlo,baz]=sphericalfwd(evla,evlo,gcarc,az)
%SPHERICALFWD    Finds a point on a sphere relative to another point
%
%    Description: [LAT2,LON2,BAZ]=SPHERICALFWD(LAT1,LON1,GCARC,AZ) returns
%     geocentric latitudes LAT2 and longitudes LON2 of destination
%     point(s), as well as the back azimuths BAZ, given the great circle
%     distances GCARC and forward azimuths AZ from initial point(s) with
%     geocentric latitudes LAT1 and longitudes LON1.  Inputs must all be in
%     degrees.  Outputs are also all in degrees.  LAT1 and LON1 must be
%     nonempty same-size arrays and GCARC and AZ must be as well.  If
%     multiple initial points and distance-azimuths are given, all must be
%     the same size (1 initial point per distance-azimuth).  A single
%     initial point may be paired with multiple distance-azimuths and
%     multiple initial points may be paired with a single distance-azimuth
%     to make working with repetitive data simpler.
%
%    Notes:
%     - Latitudes are geocentric (0 deg lat == equator, range -90<=lat<=90)
%     - Longitudes are returned in the range -180<lon<=180
%     - Azimuths are returned in the range 0<=az<=360
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    [stla,stlo]=sphericalfwd(evla,evlo,gcarc,az)
%
%    Examples:
%     St. Louis, MO USA to ???:
%      [lat,lon]=sphericalfwd(38.649,-90.305,45,-30)
%
%    See also: sphericalinv, haversine, vincentyfwd, vincentyinv

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 06:10 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin))

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(gcarc) || ~isnumeric(az))
    error('SAClab:sphericalfwd:nonNumeric','All inputs must be numeric!');
elseif(isempty(gcarc) || ~isequal(size(gcarc),size(az)))
    error('SAClab:sphericalfwd:unpairedLatLon',...
        'Distances & azimuths must be nonempty, equal size arrays!')
elseif(isempty(evla) || ~isequal(size(evla),size(evlo)))
    error('SAClab:sphericalfwd:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(~isscalar(evla) && ~isscalar(az) && ~isequal(size(evla),size(az)))
    error('SAClab:sphericalfwd:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!')
end

% convert to radians
d2r=pi/180; r2d=180/pi;
evla=evla.*d2r;
evlo=evlo.*d2r;
gcarc=gcarc.*d2r;
az=az.*d2r;

% get destination point
stla=asin(sin(evla).*cos(gcarc)+cos(evla).*sin(gcarc).*cos(az));
stlo=mod(evlo+atan2(sin(az).*sin(gcarc).*cos(evla),...
    cos(gcarc)-sin(evla).*sin(stla)),2*pi);

% get back azimuth
baz=mod(r2d.*atan2(sin(evlo-stlo).*cos(evla),...
        cos(stla).*sin(evla)-sin(stla).*cos(evla).*cos(evlo-stlo)),360);

% proper units/range
stla=stla.*r2d;
stlo=stlo.*r2d;
stlo(stlo>180)=stlo(stlo>180)-360;

end
