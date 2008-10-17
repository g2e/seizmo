function [dist,az,baz]=vincentyinv(evla,evlo,stla,stlo,ellipsoid)
%VINCENTYINV    Find distance and azimuth between 2 locations on ellipsoid
%    
%    Description: [DIST,AZ,BAZ]=VINCENTYINV(LAT1,LON1,LAT2,LON2) returns
%     the geodesic lengths DIST, forward azimuths AZ and back azimuths BAZ
%     between initial point(s) with geodetic latitudes LAT1 and longitudes
%     LON1 and final point(s) with geodetic latitudes LAT2 and longitudes
%     LON2.  All inputs must be in degrees.  DIST is in kilometers.  AZ and
%     BAZ are in degrees.  LAT1 and LON1 must be nonempty same-size arrays
%     and LAT2 and LON2 must be as well.  If multiple initial and final
%     points are given, all must be the same size (1 initial point per
%     final point).  A single initial or final point may be paired with an
%     array of the other to calculate the relative position of multiple
%     points against a single point.
%
%     VINCENTYINV(LAT1,LON1,LAT2,LON2,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  This is compatible with Matlab's Mapping Toolbox
%     function ALMANAC.
%
%    Notes:
%     - Distances and azimuths are found following the formulation of:
%        T. Vincenty (1975), Direct and Inverse Solutions of Geodesics on
%        the Ellipsoid with Application of Nested Equations, Survey Review,
%        Vol. XXII, No. 176, pp. 88-93.
%       and assume the reference ellipsoid WGS-84 unless another is given.
%     - Azimuths are returned in the range 0<=az<=360
%     - Azimuths & distances for near antipodal points are not accurate and
%       should be avoided if possible.
%
%    System requirements: Matlab 7
%
%    Usage:    [dist,az,baz]=vincentyinv(evla,evlo,stla,stlo)
%              [dist,az,baz]=vincentyinv(evla,evlo,stla,stlo,[a f])
%
%    Examples:
%     St. Louis, MO USA to Yaounde, Cameroon:
%      [dist,az,baz]=vincentyinv(38.649,-90.305,3.861,11.521)
%
%     St. Louis, USA to Isla Isabella, Galapagos:
%      [dist,az,baz]=sphericalinv(38.649,-90.305,-0.823,-91.097)
%
%    See also: vincentyfwd, sphericalinv, sphericalfwd, haversine

%     Version History:
%        June 19, 2008 - initial version
%        June 22, 2008 - handle multiple locations, specify ellipsoid
%        Oct. 11, 2008 - rename from DELAZ to VINCENTYINV, remove gcarc
%        Oct. 14, 2008 - vectorized and input size preserved
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2008 at 08:10 GMT

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
        error('SAClab:vincentyinv:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_radius flattening(<1)]'])
    end
end

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('SAClab:vincentyinv:nonNumeric','All inputs must be numeric!');
elseif(isempty(stla) || ~isequal(size(stla),size(stlo)))
    error('SAClab:vincentyinv:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(isempty(evla) || ~isequal(size(evla),size(evlo)))
    error('SAClab:vincentyinv:unpairedLatLon',...
        'Latitudes & longitudes must be nonempty, equal size arrays!')
elseif(~isscalar(evla) && ~isscalar(stla) ...
        && ~isequal(size(evla),size(stla)))
    error('SAClab:vincentyinv:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!')
end

% number of pairs
n1=size(evla); n2=size(stla); n=n1;

% expand scalars
if(isscalar(evla)); evla=repmat(evla,n2); evlo=repmat(evlo,n2); n=n2;
elseif(isscalar(stla)); stla=repmat(stla,n1); stlo=repmat(stlo,n1);
end

% check lats are within -90 to 90 (Geodetic Latitude φ)
if(any(abs(evla)>90))
    error('SAClab:vincentyinv:latitudeOutOfRange',...
        'Starting latitude out of range (-90 to 90): %g',evla);
elseif(any(abs(stla)>90))
    error('SAClab:vincentyinv:latitudeOutOfRange',...
        'Final latitude out of range (-90 to 90): %g',stla);
end

% force lon 0 to 360
evlo=mod(evlo,360);
stlo=mod(stlo,360);

% conversion constant
R2D=180/pi;
D2R=pi/180;

% reduced latitudes (in radians)
% tanU=(1-f)tanφ
U1=atan((1-f).*tan(evla.*D2R));
U2=atan((1-f).*tan(stla.*D2R));

% minimum longitude difference (in radians)
% east is positive
L=(stlo-evlo)./R2D;
L(L>pi)=L(L>pi)-2*pi;
L(L<-pi)=L(L<-pi)+2*pi;

% for efficiency
cosU1=cos(U1); cosU2=cos(U2);
sinU1=sin(U1); sinU2=sin(U2);

% various ellipsoid measures
b=a-a*f;
b2=b^2;
a2=a^2;

% iterate until λ converges
% λ is the longitude difference on an auxiliary sphere
% includes special catch for near antipodal pairs that 
% returns after first loop to avoid divergence
sigma=nan(n); sinsigma=sigma; cossigma=sigma;
sinalpha=sigma; cos2alpha=sigma; cos2sigmam=sigma; C=sigma;
lamda=L; lamdaprime=2*pi.*ones(n); eps=1e-12; oneqtr=false(n);
left=(abs(lamda-lamdaprime)>eps & abs(lamda)<(pi+eps));
while (any(left))
    % sigma = angular distance on sphere from P1 to P2
    % sigmam = angular distance on sphere from equator to midpoint
    % alpha = azimuth of the geodesic at the equator
    sinsigma(left)=sqrt((cosU2(left).*sin(lamda(left))).^2 ...
        +(cosU1(left).*sinU2(left)...
        -sinU1(left).*cosU2(left).*cos(lamda(left))).^2);
    cossigma(left)=sinU1(left).*sinU2(left)...
        +cosU1(left).*cosU2(left).*cos(lamda(left));
    sigma(left)=atan2(sinsigma(left),cossigma(left));
    sinalpha(left)=cosU1(left).*cosU2(left)...
        .*sin(lamda(left))./sinsigma(left);
    cos2alpha(left)=1-sinalpha(left).^2;
    % cos²α will be 0 for equatorial paths
    oneqtr(left)=(cos2alpha(left)==0);
    cos2sigmam(left & oneqtr)=0;
    cos2sigmam(left & ~oneqtr)=cossigma(left & ~oneqtr)...
        -2.*sinU1(left & ~oneqtr).*sinU2(left & ~oneqtr)...
        ./cos2alpha(left & ~oneqtr);
    C(left)=f/16.*cos2alpha(left).*(4+f.*(4-3.*cos2alpha(left)));
    lamdaprime(left)=lamda(left);
    lamda(left)=L(left)+(1-C(left)).*f.*sinalpha(left).*(sigma(left)...
        +C(left).*sinsigma(left).*(cos2sigmam(left)...
        +C(left).*cossigma(left).*(-1+2.*cos2sigmam(left).^2)));
    left(left)=(abs(lamda(left)-lamdaprime(left))>eps ...
        & abs(lamda(left))<(pi+eps));
end

% get the distance and azimuth
u2=cos2alpha.*(a2-b2)./b2;
A=1+u2./16384.*(4096+u2.*(-768+u2.*(320-175.*u2)));
B=u2./1024.*(256+u2.*(-128+u2.*(74-47.*u2)));
deltasigma=B.*sinsigma.*(cos2sigmam...
    +B./4.*(cossigma.*(-1+2.*cos2sigmam.^2)...
    -B./6.*cos2sigmam.*(-3+4.*sinsigma.^2).*(-3+4.*cos2sigmam.^2)));
dist=b.*A.*(sigma-deltasigma)/1000;
az=mod(atan2(cosU2.*sin(lamda),...
    cosU1.*sinU2-sinU1.*cosU2.*cos(lamda)).*R2D,360);
baz=mod(atan2(-cosU1.*sin(lamda),...
    sinU1.*cosU2-cosU1.*sinU2.*cos(lamda)).*R2D,360);

end
