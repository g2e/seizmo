function [gcarc]=haversine(evla,evlo,stla,stlo)
%HAVERSINE    Returns distance between 2 points using the Haversine formula
%
%    Description: HAVERSINE(LAT1,LON1,LAT2,LON2) returns the spherical
%     greater-circle-arc degree distance between two points.  LAT1, LON1,
%     LAT2, LON2 must all be in degrees and the latitudes must be
%     geocentric.  LAT1 and LON1 must be nonempty same-size arrays and LAT2
%     and LON2 must be as well.  If multiple initial and final points are
%     given, they must be the same size (1 initial point per final point).
%     A single initial or final point may be paired with an array of the
%     other to calculate the distance of multiple points from a single
%     point.
%
%    Notes:
%     - 'half versed sine' is better suited for accuracy at small distances
%       than SPHERICALINV as it uses the haversine function rather than a
%       cosine which becomes inefficient at small distances.
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    gcarc=haversine(lat1,lon1,lat2,lon2)
%
%    Examples:
%
%    See also: sphericalinv, vincentyinv, sphericalfwd, vincentyfwd

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 12, 2008 at 17:00 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin))

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('SAClab:haversine:nonNumeric','All inputs must be numeric!');
elseif(isempty(stla) || ~isequal(size(stla),size(stlo)))
    error('SAClab:haversine:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(isempty(evla) || ~isequal(size(evla),size(evlo)))
    error('SAClab:haversine:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(~isscalar(evla) && ~isscalar(stla) ...
        && ~isequal(size(evla),size(stla)))
    error('SAClab:haversine:nonscalarUnequalArrays',...
        'Location arrays need to be scalar or have equal size!')
end

% get haversine distance
a=sind((stla-evla)/2).^2+cosd(evla).*cosd(stla).*sind((stlo-evlo)/2).^2;
gcarc=2.*atan2(sqrt(a),sqrt(1-a)).*(180/pi);

end
