function [stla,stlo,baz]=sphericalfwd(evla,evlo,gcarc,az)
%SPHERICALFWD    Finds a point on a sphere relative to another point
%
%    Usage:    [lat2,lon2,baz]=sphericalfwd(lat1,lon1,gcarc,az)
%
%    Description:
%     [LAT2,LON2,BAZ]=SPHERICALFWD(LAT1,LON1,GCARC,AZ) returns geocentric
%     latitudes LAT2 and longitudes LON2 of destination point(s), as well
%     as the back azimuths BAZ, given the great circle distances GCARC and
%     forward azimuths AZ from initial point(s) with geocentric latitudes
%     LAT1 and longitudes LON1.  Inputs must all be in degrees.  Outputs
%     are also all in degrees.
%
%    Notes:
%     - Latitudes are geocentric (0 deg lat == equator, range +/-90)
%     - Longitudes are returned in the range +/-180
%     - Backazimuth is returned in the range 0-360
%
%    Examples:
%     % St. Louis, MO USA to ???:
%     [lat,lon,baz]=sphericalfwd(38.649,-90.305,45,-30)
%
%    See also: SPHERICALINV, HAVERSINE, VINCENTYFWD, VINCENTYINV

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - improved scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        May   8, 2009 - minor doc fix
%        Jan. 22, 2011 - use degrees functions, fix pole result giving
%                        complex occasionally, nargchk fix
%        Feb. 24, 2011 - doc update
%        Feb.  9, 2012 - doc update, drop replication for speed, precompute
%                        sines for speed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 21:15 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin));

% size up inputs
sz{1}=size(evla); sz{2}=size(evlo);
sz{3}=size(gcarc); sz{4}=size(az);
n(1)=prod(sz{1}); n(2)=prod(sz{2});
n(3)=prod(sz{3}); n(4)=prod(sz{4});

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(gcarc) || ~isnumeric(az))
    error('seizmo:sphericalfwd:nonNumeric','All inputs must be numeric!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:sphericalfwd:badSize',...
        'All inputs must be equal sized or scalar!');
end

% optimize computation over memory
sinlat1=sind(evla);
singc=sind(gcarc);
coslat1=cosd(evla);
cosgc=cosd(gcarc);

% for conversion
r2d=180/pi;

% get destination point
stla=real(asind(sinlat1.*cosgc+coslat1.*singc.*cosd(az)));
sinlat2=sind(stla);
stlo=mod(evlo+r2d.*atan2(sind(az).*singc.*coslat1,...
    cosgc-sinlat1.*sinlat2),360);

% stla misses evlo for sizing
if(isscalar(stla)); stla=stla(ones(sz{2})); end

% proper range (+/-180)
stlo(stlo>180)=stlo(stlo>180)-360;

% shortcut
if(nargout<3); return; end

% get back azimuth
baz=mod(r2d.*atan2(sind(evlo-stlo).*coslat1,...
        cosd(stla).*sinlat1-sinlat2.*coslat1.*cosd(evlo-stlo)),360);

end
