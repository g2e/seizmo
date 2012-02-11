function [gcarc]=haversine(evla,evlo,stla,stlo)
%HAVERSINE    Returns distance between 2 points using the Haversine formula
%
%    Usage:    gcarc=haversine(lat1,lon1,lat2,lon2)
%
%    Description:
%     HAVERSINE(LAT1,LON1,LAT2,LON2) returns the spherical greater-circle-
%     arc degree distance between two points.  LAT1, LON1, LAT2, LON2 must
%     all be in degrees and the latitudes must be geocentric.
%
%    Notes:
%     - 'half versed sine' is better suited for accuracy at small distances
%       compared to SPHERICALINV as it uses the haversine function rather
%       than a cosine which becomes inefficient at small distances.
%
%    Examples:
%     % Plot distance discrepancy for sphericalinv and haversine:
%     deg2m=1000*6371*pi/180;
%     dist=10.^(-1:.01:7)'./deg2m;
%     loglog(deg2m*dist,deg2m*abs(dist-sphericalinv(0,0,dist,0)),...
%            deg2m*dist,deg2m*abs(dist-haversine(0,0,dist,0)))
%     ylabel('discrepancy (m)')
%     xlabel('distance (m)')
%     legend({'sphericalinv' 'haversine'})
%     % demonstrates convincingly this function is more accurate!
%
%    See also: SPHERICALINV, VINCENTYINV, SPHERICALFWD, VINCENTYFWD

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - improved scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  9, 2012 - doc update, drop replication for speed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 15:05 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin));

% size up inputs
sz{1}=size(evla); sz{2}=size(evlo);
sz{3}=size(stla); sz{4}=size(stlo);
n(1)=prod(sz{1}); n(2)=prod(sz{2});
n(3)=prod(sz{3}); n(4)=prod(sz{4});

% check inputs
if(~isreal(evla) || ~isreal(evlo) ||...
        ~isreal(stla) || ~isreal(stlo))
    error('seizmo:haversine:nonNumeric',...
        'All inputs must be real-valued arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:haversine:badSize',...
        'All inputs must be equal sized or scalar!');
end

% get haversine distance
a=sind((stla-evla)/2).^2+cosd(evla).*cosd(stla).*sind((stlo-evlo)/2).^2;
gcarc=2.*atan2(sqrt(a),sqrt(1-a)).*(180/pi);

end
