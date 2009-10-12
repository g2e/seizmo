function [gcarc,az,baz]=sphericalinv(evla,evlo,stla,stlo)
%SPHERICALINV    Return distance and azimuth between 2 locations on sphere
%
%    Usage:    [gcarc,az,baz]=sphericalinv(lat1,lon1,lat2,lon2)
%
%    Description: [GCARC,AZ,BAZ]=SPHERICALINV(LAT1,LON1,LAT2,LON2) returns
%     the great-circle-arc degree distances GCARC, forward azimuths AZ and
%     back azimuths BAZ between initial point(s) with geocentric latitudes
%     LAT1 and longitudes LON1 and final point(s) with geocentric latitudes
%     LAT2 and longitudes LON2 on a sphere.  All inputs must be in degrees.
%     Outputs are also all in degrees.  LAT1 and LON1 must be scalar or
%     nonempty same-size arrays and LAT2 and LON2 must be as well.  If
%     multiple initial and final points are given, all must be the same
%     size (1 initial point per final point).  A single initial or final
%     point may be paired with an array of the other to calculate the
%     relative position of multiple points against a single point.
%
%    Notes:
%     - Will always return the shorter great-circle-arc (GCARC<=180)
%     - Accuracy degrades at very small distances (see HAVERSINE)
%     - Azimuths are returned in the range 0<=az<=360
%
%    Examples:
%     St. Louis, MO USA to Yaounde, Cameroon:
%      [dist,az,baz]=sphericalinv(38.649,-90.305,3.861,11.521)
%
%     St. Louis, MO USA to Isla Isabella, Galapagos:
%      [dist,az,baz]=sphericalinv(38.649,-90.305,-0.823,-91.097)
%
%    See also: HAVERSINE, SPHERICALFWD, VINCENTYINV, VINCENTYFWD

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - improved scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% require 4 inputs
msg=nargchk(4,4,nargin);
if(~isempty(msg)); error(msg); end

% size up inputs
sz1=size(evla); sz2=size(evlo);
sz3=size(stla); sz4=size(stlo);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('seizmo:sphericalinv:nonNumeric','All inputs must be numeric!');
elseif(any([n1 n2 n3 n4]==0))
    error('seizmo:sphericalinv:emptyLatLon',...
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
    error('seizmo:sphericalinv:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n2==1); evla=repmat(evla,sz3); evlo=repmat(evlo,sz3); end
if(n4==1); stla=repmat(stla,sz1); stlo=repmat(stlo,sz1); end

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
