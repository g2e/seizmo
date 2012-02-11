function [dist,az,baz]=vincentyinv(evla,evlo,stla,stlo,ellipsoid,tolerance)
%VINCENTYINV    Find distance and azimuth between 2 locations on ellipsoid
%    
%    Usage:    [dist,az,baz]=vincentyinv(lat1,lon1,lat2,lon2)
%              [dist,az,baz]=vincentyinv(lat1,lon1,lat2,lon2,[a f])
%              [dist,az,baz]=vincentyinv(lat1,lon1,lat2,lon2,[a f],...
%                                        [tol1 tol2])
%
%    Description:
%     [DIST,AZ,BAZ]=VINCENTYINV(LAT1,LON1,LAT2,LON2) returns the geodesic
%     lengths DIST, forward azimuths AZ and back azimuths BAZ between
%     initial point(s) with geographic latitudes LAT1 & longitudes LON1 and
%     final point(s) with geographic latitudes LAT2 & longitudes LON2 on
%     the WGS-84 reference ellipsoid.  All inputs must be in degrees.  DIST
%     is in kilometers.  AZ and BAZ are in degrees.
%
%     VINCENTYINV(LAT1,LON1,LAT2,LON2,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  This is compatible with output from Matlab's Mapping
%     Toolbox function ALMANAC.  By default the ellipsoid parameters are
%     set to those of the WGS-84 reference ellipsoid.
%
%     VINCENTYINV(LAT1,LON1,LAT2,LON2,[A F],[TOL1 TOL2]) allows specifying
%     the tolerances of the calculation.  TOL1 sets the minimum precision
%     in radians for convergent cases.  TOL2 sets the maximum distance in
%     radians from the antipode at which a case is deemed divergent and is
%     forced to return with the latest result.  Note that setting TOL2 too
%     high can force convergent cases to return before they satify TOL1.
%     By default T0L1 is 1e-12 and TOL2 is 1e-2.
%
%    Notes:
%     - Distances and azimuths are found following the formulation of:
%        T. Vincenty (1975), Direct and Inverse Solutions of Geodesics on
%        the Ellipsoid with Application of Nested Equations, Survey Review,
%        Vol. XXII, No. 176, pp. 88-93.
%       and assume the reference ellipsoid WGS-84 unless another is given.
%     - Azimuths are returned in the range 0-360
%     - Distances and, in particular, azimuths for near antipodal points
%       are not accurate and should be avoided if possible.
%
%    Examples:
%     % St. Louis, MO USA to Yaounde, Cameroon:
%     [dist,az,baz]=vincentyinv(38.649,-90.305,3.861,11.521)
%
%     % St. Louis, MO USA to Isla Isabella, Galapagos:
%     [dist,az,baz]=sphericalinv(38.649,-90.305,-0.823,-91.097)
%
%    See also: VINCENTYFWD, SPHERICALINV, SPHERICALFWD, HAVERSINE

%     Version History:
%        June 19, 2008 - initial version
%        June 22, 2008 - handle multiple locations, specify ellipsoid
%        Oct. 11, 2008 - rename from DELAZ to VINCENTYINV, remove gcarc
%        Oct. 14, 2008 - vectorized and input size preserved
%        Oct. 26, 2008 - fix infinite loop for near antipode cases, allow
%                        specifying the tolerance for convergent and
%                        divergent cases, improved scalar expansion, doc
%                        and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        May   8, 2009 - minor doc fix
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Apr. 10, 2010 - return 0s rather than NaNs for equal positions
%        Feb. 10, 2011 - minor code formating changes
%        Feb.  9, 2012 - doc update, drop some replication for speed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 13:15 GMT

% todo:

% require 4 to 6 inputs
error(nargchk(4,6,nargin));

% default ellipsoid - WGS-84 Reference Ellipsoid
if(nargin==4 || isempty(ellipsoid))
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
        error('seizmo:vincentyinv:badEllipsoid',...
            ['Ellipsoid must be a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% default tolerances
if(any(nargin==[4 5]) || isempty(tolerance))
    % ~0.01mm and ~100km
    tolerance(1:2)=[1e-12 1e-2];
elseif(~isnumeric(tolerance) || numel(tolerance)>2)
    error('seizmo:vincentyinv:badTolerance',...
        ['Tolerance must be a 2 element vector specifying:\n'...
        '[precision antipodal_precision]']);
elseif(any(tolerance<0) || any(tolerance>pi))
    error('seizmo:vincentyinv:badTolerance',...
        'Tolerances must be in the range 0 to PI!');
elseif(numel(tolerance)==1)
    tolerance(2)=tolerance(1);
end

% size up inputs
sz{1}=size(evla); sz{2}=size(evlo);
sz{3}=size(stla); sz{4}=size(stlo);
n(1)=prod(sz{1}); n(2)=prod(sz{2});
n(3)=prod(sz{3}); n(4)=prod(sz{4});

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('seizmo:vincentyinv:nonNumeric','All inputs must be numeric!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:vincentyinv:badSize',...
        'All location inputs must be equal sized or scalar!');
end

% equal location check (also expanded)
eqpo=(evla==stla & evlo==stlo);

% check lats are within -90 to 90 (Geographic Latitude φ)
if(any(abs(evla)>90))
    error('seizmo:vincentyinv:latitudeOutOfRange',...
        'Starting latitudes out of range (-90 to 90)');
elseif(any(abs(stla)>90))
    error('seizmo:vincentyinv:latitudeOutOfRange',...
        'Final latitudes out of range (-90 to 90)');
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

% pole latitude measure
% 0 for antipodal or equator pairs
% PI^2 for pairs at a pole
PLM=abs(U1+U2).^2;

% minimum longitude difference (in radians)
% east is positive
L=(stlo-evlo).*D2R;
L(L>pi)=L(L>pi)-2*pi;
L(L<-pi)=L(L<-pi)+2*pi;

% for efficiency
cosU1=cos(U1); cosU2=cos(U2);
sinU1=sin(U1); sinU2=sin(U2);

% scalar expansion (needed for loop)
n=size(eqpo);
if(~isequal(sz{:}))
    if(isscalar(PLM)); PLM=PLM(ones(n)); end
    if(isscalar(L)); L=L(ones(n)); end
    if(isscalar(cosU1)); cosU1=cosU1(ones(n)); end
    if(isscalar(cosU2)); cosU2=cosU2(ones(n)); end
    if(isscalar(sinU1)); sinU1=sinU1(ones(n)); end
    if(isscalar(sinU2)); sinU2=sinU2(ones(n)); end
end

% various ellipsoid measures
b=a-a*f;
b2=b^2;
a2=a^2;

% iterate until λ converges or is deemed divergent
% λ is the longitude difference on an auxiliary sphere
% includes special catch for near antipodal pairs
sigma=nan(n); sinsigma=sigma; cossigma=sigma;
sinalpha=sigma; cos2alpha=sigma; cos2sigmam=sigma; C=sigma;
lamda=L; lamdaprime=2*pi.*ones(n); on_eqtr=false(n); left=true(n);
while (any(left)) % forced through at least one iteration
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
    on_eqtr(left)=(cos2alpha(left)==0);
    cos2sigmam(left & on_eqtr)=0;
    cos2sigmam(left & ~on_eqtr)=cossigma(left & ~on_eqtr)...
        -2.*sinU1(left & ~on_eqtr).*sinU2(left & ~on_eqtr)...
        ./cos2alpha(left & ~on_eqtr);
    C(left)=f/16.*cos2alpha(left).*(4+f.*(4-3.*cos2alpha(left)));
    lamdaprime(left)=lamda(left);
    lamda(left)=L(left)+(1-C(left)).*f.*sinalpha(left).*(sigma(left)...
        +C(left).*sinsigma(left).*(cos2sigmam(left)...
        +C(left).*cossigma(left).*(-1+2.*cos2sigmam(left).^2)));
    left(left)=(abs(lamda(left)-lamdaprime(left))>tolerance(1)...
        & abs(lamda(left))<=pi...
        & ~(sqrt((pi-abs(lamda(left))).^2+PLM(left))<tolerance(2)));
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

% avoids NaNs when getting geometry between 2 equal positions
if(any(eqpo(:))); dist(eqpo)=0; az(eqpo)=0; baz(eqpo)=0; end

end
