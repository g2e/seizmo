function [gcarc]=haversine(evla,evlo,stla,stlo)
%HAVERSINE    Returns distance between 2 points using the Haversine formula
%
%    Usage:    gcarc=haversine(lat1,lon1,lat2,lon2)
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
%    Examples:
%     Plotting up short distance results for sphericalinv and haversine:
%      plot(0:1e-9:1e-5,sphericalinv(0,0,0:1e-9:1e-5,0),...
%           0:1e-9:1e-5,haversine(0,0,0:1e-9:1e-5,0))
%     demonstrates where this function becomes useful (couple meters).
%
%    See also: SPHERICALINV, VINCENTYINV, SPHERICALFWD, VINCENTYFWD

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - improved scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin));

% size up inputs
sz1=size(evla); sz2=size(evlo);
sz3=size(stla); sz4=size(stlo);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('seizmo:haversine:nonNumeric','All inputs must be numeric!');
elseif(any([n1 n2 n3 n4]==0))
    error('seizmo:haversine:emptyLatLon',...
        'Latitudes & longitudes must be nonempty arrays!');
end

% expand scalars
if(n1==1); evla=repmat(evla,sz2); n1=n2; sz1=sz2; end
if(n2==1); evlo=repmat(evlo,sz1); n2=n1; sz2=sz1; end
if(n3==1); stla=repmat(stla,sz4); n3=n4; sz3=sz4; end
if(n4==1); stlo=repmat(stlo,sz3); n4=n3; sz4=sz3; end

% cross check inputs
if(~isequal(sz1,sz2) || ~isequal(sz3,sz4) ||...
        (~any([n1 n3]==1) && ~isequal(sz1,sz3)))
    error('seizmo:haversine:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n2==1); evla=repmat(evla,sz3); evlo=repmat(evlo,sz3); end
if(n4==1); stla=repmat(stla,sz1); stlo=repmat(stlo,sz1); end

% get haversine distance
a=sind((stla-evla)/2).^2+cosd(evla).*cosd(stla).*sind((stlo-evlo)/2).^2;
gcarc=2.*atan2(sqrt(a),sqrt(1-a)).*(180/pi);

end
