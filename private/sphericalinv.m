function [gcarc,az,baz]=sphericalinv(evla,evlo,stla,stlo)
%SPHERICALINV    Returns great-circle dist,azim between 2 points on sphere
%
%    Description: [GCARC,AZ,BAZ]=SPHERICALINV(LAT1,LON1,LAT2,LON2) returns
%     the great-circle-arc degree distances GCARC, forward azimuths AZ and
%     back azimuths BAZ between initial point(s) with geocentric latitudes
%     LAT1 and longitudes LON1 and final point(s) with geocentric latitudes
%     LAT2 and longitudes LON2.  All inputs must be in degrees.  Outputs
%     are also all in degrees.  LAT1 and LON1 must be nonempty same-size
%     arrays and LAT2 and LON2 must be as well.  If multiple initial and
%     final points are given, all must be the same size (1 initial point
%     per final point).  A single initial or final point may be paired with
%     an array of the other to calculate the relative position of multiple
%     points against a single point.
%
%    Notes:
%     - Will always return the shorter great-circle-arc (GCARC<=180)
%     - Accuracy degrades at very small distances (see HAVERSINE)
%     - Azimuths are returned in the range 0<=az<=360
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    [gcarc,az,baz]=sphericalinv(lat1,lon1,lat2,lon2)
%
%    Examples:
%     St. Louis, USA to Yaounde, Cameroon:
%      [dist,az,baz]=sphericalinv(38.649,-90.305,3.861,11.521)
%
%     St. Louis, USA to Isla Isabella, Galapagos:
%      [dist,az,baz]=sphericalinv(38.649,-90.305,-0.823,-91.097)
%
%    See also: haversine, sphericalfwd, vincentyinv, vincentyfwd

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 06:45 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin))

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('SAClab:sphericalinv:nonNumeric','All inputs must be numeric!');
elseif(isempty(stla) || ~isequal(size(stla),size(stlo)))
    error('SAClab:sphericalinv:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(isempty(evla) || ~isequal(size(evla),size(evlo)))
    error('SAClab:sphericalinv:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(~isscalar(evla) && ~isscalar(stla) ...
        && ~isequal(size(evla),size(stla)))
    error('SAClab:sphericalinv:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!')
end

% convert to radians
d2r=pi/180; r2d=1/d2r;
evla=evla.*d2r;
evlo=evlo.*d2r;
stla=stla.*d2r;
stlo=stlo.*d2r;

% optimize
sinlat1=sin(evla);
sinlat2=sin(stla);
coslat1=cos(evla);
coslat2=cos(stla);
coslo=cos(stlo-evlo);
sinlo=sin(stlo-evlo);

% get law-of-cosines distance
gcarc=acosd(sinlat1.*sinlat2+coslat1.*coslat2.*coslo);

% azimuths
az=mod(r2d.*atan2(sinlo.*coslat2,...
    coslat1.*sinlat2-sinlat1.*coslat2.*coslo),360);
baz=mod(r2d.*atan2(-sinlo.*coslat1,...
    coslat2.*sinlat1-sinlat2.*coslat1.*coslo),360);

end
