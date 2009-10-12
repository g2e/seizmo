function [stla,stlo,baz]=sphericalfwd(evla,evlo,gcarc,az)
%SPHERICALFWD    Finds a point on a sphere relative to another point
%
%    Usage:    [lat2,lon2]=sphericalfwd(lat1,lon1,gcarc,az)
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
%    Examples:
%     St. Louis, MO USA to ???:
%      [lat,lon]=sphericalfwd(38.649,-90.305,45,-30)
%
%    See also: SPHERICALINV, HAVERSINE, VINCENTYFWD, VINCENTYINV

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - improved scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        May   8, 2009 - minor doc fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% require 4 inputs
msg=nargchk(4,4,nargin);
if(~isempty(msg)); error(msg); end

% size up inputs
sz1=size(evla); sz2=size(evlo);
sz3=size(gcarc); sz4=size(az);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(gcarc) || ~isnumeric(az))
    error('seizmo:sphericalfwd:nonNumeric','All inputs must be numeric!');
elseif(any([n1 n2 n3 n4]==0))
    error('seizmo:sphericalfwd:emptyLatLon',...
        'Location inputs must be nonempty arrays!');
end

% expand scalars
if(n1==1); evla=repmat(evla,sz2); n1=n2; sz1=sz2; end
if(n2==1); evlo=repmat(evlo,sz1); n2=n1; sz2=sz1; end
if(n3==1); stla=repmat(gcarc,sz4); n3=n4; sz3=sz4; end
if(n4==1); stlo=repmat(az,sz3); n4=n3; sz4=sz3; end

% cross check inputs
if(~isequal(sz1,sz2) || ~isequal(sz3,sz4) ||...
        (~any([n1 n3]==1) && ~isequal(sz1,sz3)))
    error('seizmo:sphericalfwd:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n2==1); evla=repmat(evla,sz3); evlo=repmat(evlo,sz3); end
if(n4==1); gcarc=repmat(stla,sz1); az=repmat(stlo,sz1); end

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
