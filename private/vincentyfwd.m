function [stla,stlo,baz]=vincentyfwd(evla,evlo,dist,az,ellipsoid)
%VINCENTYFWD    Find destination point on an ellipsoid relative to a point
%    
%    Description: [LAT2,LON2,BAZ]=VINCENTYFWD(LAT1,LON1,DIST,AZ) returns
%     geodetic latitudes LAT2 and longitudes LON2 of destination point(s),
%     as well as the backazimuths BAZ, given the distances DIST and forward
%     azimuths AZ from initial point(s) with geodetic latitudes LAT1 and
%     longitudes LON1.  Inputs are all in degrees except DIST which must be
%     in kilometers.  Outputs are all in degrees.  LAT1 and LON1 must be
%     nonempty same-size arrays and DIST and AZ must be as well.  If
%     multiple initial points and distance-azimuths are given, all must be
%     same size (1 initial point per distance-azimuth).  A single initial
%     point may be paired with multiple distance-azimuths and multiple
%     initial points may be paired with a single distance-azimuth to make
%     working with repetitive data simpler.
%
%     VINCENTYFWD(LAT1,LON1,DIST,AZ,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  This is compatible with Matlab's Mapping Toolbox
%     function ALMANAC.
%
%    Notes:
%     - Destination points are found following the formulation of:
%        T. Vincenty (1975), Direct and Inverse Solutions of Geodesics on
%        the Ellipsoid with Application of Nested Equations, Survey Review,
%        Vol. XXII, No. 176, pp. 88-93.
%       and assume the reference ellipsoid WGS-84 unless another is given.
%     - Latitudes are geodetic (0 deg lat == equator, range -90<=lat<=90)
%     - Longitudes are returned in the range -180<lon<=180
%     - Azimuths are returned in the range 0<=az<=360
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    [lat2,lon2,baz]=vincentyfwd(lat1,lon1,dist,az)
%              [lat2,lon2,baz]=vincentyfwd(lat1,lon1,dist,az,[a f])
%
%    Examples:
%     St. Louis, MO USA to ???:
%      [lat2,lon2,baz]=vincentyfwd(38.649,-90.305,5000,-30)
%
%    See also: vincentyinv, sphericalinv, sphericalfwd, haversine

%     Version History:
%        Oct. 14, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 05:45 GMT

% todo:

% require 4 or 5 inputs
error(nargchk(4,5,nargin))

% default - WGS-84 Reference Ellipsoid
if(nargin==4)
    % a=radius at equator (major axis)
    % f=flattening
    a=6378137.0;
    f=1/298.257223563;
else
    % manually specify ellipsoid (will accept almanac output)
    if(isnumeric(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1)*1000; % km=>m
        f=ellipsoid(2);
    else
        error('SAClab:vincentyfwd:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_radius flattening(<1)]'])
    end
end

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(dist) || ~isnumeric(az))
    error('SAClab:vincentyfwd:nonNumeric','All inputs must be numeric!');
elseif(isempty(dist) || ~isequal(size(dist),size(az)))
    error('SAClab:vincentyfwd:unpairedLatLon',...
        'Distances and azimuths must be nonempty, equal size arrays!')
elseif(isempty(evla) || ~isequal(size(evla),size(evlo)))
    error('SAClab:vincentyfwd:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(~isscalar(evla) && ~isscalar(az) ...
        && ~isequal(size(evla),size(az)))
    error('SAClab:vincentyfwd:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!')
end

% number of pairs
n1=size(evla); n2=size(dist); n=n1;

% expand scalars
if(isscalar(evla)); evla=repmat(evla,n2); evlo=repmat(evlo,n2); n=n2;
elseif(isscalar(az)); dist=repmat(dist,n1); az=repmat(az,n1);
end

% check lats are within -90 to 90 (Geodetic Latitude φ)
if(any(abs(evla)>90))
    error('SAClab:vincentyfwd:latitudeOutOfRange',...
        'Starting latitude out of range (-90 to 90): %g',evla);
end

% force lon 0 to 360
evlo=mod(evlo,360);

% conversion constants
R2D=180/pi;
D2R=pi/180;

% convert to radians
az=az.*D2R;

% convert to meters
dist=dist.*1000;

% reduced latitude
tanU1=(1-f).*tan(evla.*D2R);
cosU1=1./sqrt(1+tanU1.^2);
sinU1=tanU1.*cosU1;

% various ellipsoid measures
b=a-a*f;
b2=b^2;
a2=a^2;

% setup
sinalpha1=sin(az); cosalpha1=cos(az);
sigma1=atan2(tanU1,cosalpha1);
sinalpha=cosU1.*sinalpha1;
sin2alpha=sinalpha.^2;
cos2alpha=1-sin2alpha;
u2=cos2alpha.*(a2-b2)./b2;
A=1+u2./16384.*(4096+u2.*(-768+u2.*(320-175.*u2)));
B=u2./1024.*(256+u2.*(-128+u2.*(74-47.*u2)));

% iterate until sigma converges
left=true(n); cos2sigmam=nan(n); deltasigma=nan(n);
sigma=dist./(b.*A); sigmaprime=2*pi*ones(n); eps=1e-12;
while (any(left)) % forces at least one iteration
    cos2sigmam(left)=cos(2.*sigma1(left)+sigma(left));
    deltasigma(left)=B(left).*sin(sigma(left)).*(cos2sigmam(left)...
        +B(left)./4.*(cos(sigma(left)).*(-1+2.*cos2sigmam(left).^2)...
        -B(left)./6.*cos2sigmam(left).*(-3+4.*sin(sigma(left)).^2)...
        .*(-3+4.*cos2sigmam(left).^2)));
    sigmaprime(left)=sigma(left);
    sigma(left)=dist(left)./(b.*A(left))+deltasigma(left);
    left(left)=abs(sigma(left)-sigmaprime(left))>eps;
end

% get destination point
cossigma=cos(sigma); sinsigma=sin(sigma); cos2sigmam=cos(2.*sigma1+sigma);
stla=R2D.*atan2(sinU1.*cossigma+cosU1.*sinsigma.*cosalpha1,...
    (1-f).*sqrt(sin2alpha+(sinU1.*sinsigma...
    -cosU1.*cossigma.*cosalpha1).^2));
lambda=atan2(sinsigma.*sinalpha1,...
    cosU1.*cossigma-sinU1.*sinsigma.*cosalpha1);
C=f./16.*cos2alpha.*(4+f.*(4-3.*cos2alpha));
L=lambda-(1-C).*f.*sinalpha.*(sigma+C.*sinsigma.*(cos2sigmam+...
    C.*cossigma.*(-1+2.*cos2sigmam.^2)));
stlo=mod(evlo+L*R2D,360);
stlo(stlo>180)=stlo(stlo>180)-360;

% get backazimuth
baz=mod(180+...
    R2D.*atan2(sinalpha,-sinU1.*sinsigma+cosU1.*cossigma.*cosalpha1),360);

end
