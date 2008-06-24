function [gcarc,az,baz,dist]=delaz(evla,evlo,stla,stlo,ellipsoid)
%DELAZ    Find distance and azimuth between two locations on ellipsoid
%    
%    Description: [GCARC,AZ,BAZ,DIST]=DELAZ(EVLA,EVLO,STLA,STLO) computes
%     the distances (angular degree distance GCARC and kilometer distance 
%     DIST) and azimuths (forward azimuth AZ and backazimuth BAZ in 
%     degrees) from one location (geodetic latitude EVLA and longitude 
%     EVLO) to another location (geodetic latitude STLA and longitude 
%     STLO).  Multiple starting locations to a single end location as well
%     as a single starting location to multiple end locations are
%     acceptable by inputting arrays for the location info.  Otherwise the
%     number of starting and ending locations must match.
%
%     DELAZ(EVLA,EVLO,STLA,STLO,[A F]) allows specifying the ellipsoid
%     parameters A (equatorial radius in kilometers) and F (flattening).
%     This allows combination with Matlab's Mapping Toolbox function
%     ALMANAC.
%
%    Notes:
%     - distances and azimuths are found following the formulation of:
%
%        T. Vincenty (1975), Direct and Inverse Solutions of Geodesics on
%        the Ellipsoid with Application of Nested Equations, Survey Review,
%        Vol. XXII, No. 176, pp. 88-93.
%
%       and assume the reference ellipsoid WGS-84.
%     - azimuths/distances for near antipodal points are not accurate
%
%    System requirements: Matlab 7
%
%    Data requirements: 
%     - Latitudes and Longitudes in Degrees 
%     - Geodetic Latitudes (-90 to 90)
%
%    Usage: [gcarc,az,baz,dist]=delaz(evla,evlo,stla,stlo)
%
%    Examples:
%     St. Louis, USA to Yaounde, Cameroon:
%      [gcarc,az,baz,dist]=delaz(38.649,-90.305,3.861,11.521)
%
%    See also:

%     Version History:
%        June 19, 2008 - Initial Version
%        June 22, 2008 - Handle multiple locations, specify ellipsoid
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 22, 2008 at 19:30 GMT

% todo:
% gcarc - is it spherical or what

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
        error('SAClab:delaz:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_radius flattening(<1)]'])
    end
end

% other useful term
% b=radius at pole (minor axis)
b=a-a*f;

% check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) || ~isnumeric(stla) || ~isnumeric(stlo))
    error('SAClab:delaz:nonNumeric','All inputs must be numeric!');
elseif(isempty(stla) || isempty(stlo) || numel(stla)~=numel(stlo))
    error('SAClab:delaz:unpairedLatLon',...
        ['Station latitudes & longitudes must be'...
        ' nonempty, equal size arrays!'])
elseif(isempty(evla) || isempty(evlo) || numel(evla)~=numel(evlo))
    error('SAClab:delaz:unpairedLatLon',...
        ['Event latitudes & longitudes must be'...
        ' nonempty, equal size arrays!'])
elseif(~isscalar(evla) && ~isscalar(stla) && numel(evla)~=numel(stla))
    error('SAClab:delaz:nonscalarUnequalArrays',...
        'Event/Station locations arrays need to be scalar or have equal size')
end

% number of event/station pairs
n=max([numel(stla) numel(evla)]);

% expand scalar locations
if(isscalar(evla)); evla(1:n,1)=evla; evlo(1:n,1)=evlo;
elseif(isscalar(stla)); stla(1:n,1)=stla; stlo(1:n,1)=stlo;
end

% check lats are within -90 to 90 (Geodetic Latitude φ)
if(any(abs(evla)>90))
    error('SAClab:delaz:latitudeOutOfRange',...
        'Event latitude out of range (-90 to 90): %g',evla);
elseif(any(abs(stla)>90))
    error('SAClab:delaz:latitudeOutOfRange',...
        'Station latitude out of range (-90 to 90): %g',stla);
end

% force lon 0 to 360
evlo=mod(evlo,360);
stlo=mod(stlo,360);

% conversion constant
R2D=180/pi;

% reduced latitudes (in radians)
% tanU=(1-f)tanφ
U1=atan((1-f)*tan(evla/R2D));
U2=atan((1-f)*tan(stla/R2D));

% minimum longitude difference (in radians)
% east is positive
L=(stlo-evlo)/R2D;
L(L>pi)=L(L>pi)-2*pi;
L(L<-pi)=L(L<-pi)+2*pi;

% for efficiency
cosU1=cos(U1); cosU2=cos(U2);
sinU1=sin(U1); sinU2=sin(U2);

% loop through each event/station pair separately
az=zeros(n,1); baz=az; gcarc=az; dist=az;
for i=1:n
    % iterate until λ converges
    % λ is the longitude difference on an auxiliary sphere
    % special catch for near antipodal pairs returns after first loop
    lamda=L(i); lamda2=2*pi; eps=1e-12;
    while (abs(lamda-lamda2)>eps && abs(lamda)<(pi+eps))
        % s = angular distance on sphere from P1 to P2
        % sm = angular distance on sphere from equator to midpoint
        % a = azimuth of the geodesic at the equator
        sins=sqrt((cosU2(i)*sin(lamda))^2+(cosU1(i)*sinU2(i)-sinU1(i)*cosU2(i)*cos(lamda))^2);
        coss=sinU1(i)*sinU2(i)+cosU1(i)*cosU2(i)*cos(lamda);
        s=atan2(sins,coss);
        sina=cosU1(i)*cosU2(i)*sin(lamda)/sins;
        cosa2=1-sina^2;
        % cos²α will be 0 for equatorial paths
        if(cosa2); cos2sm=coss-2*sinU1(i)*sinU2(i)/cosa2;
        else cos2sm=0;
        end
        C=f/16*cosa2*(4+f*(4-3*cosa2));
        lamda2=lamda;
        lamda=L(i)+(1-C)*f*sina*(s+C*sins*(cos2sm+C*coss*(-1+2*cos2sm^2)));
    end
    
    % get the distance and azimuth
    u2=cosa2*(a^2-b^2)/b^2;
    A=1+u2/16384*(4096+u2*(-768+u2*(320-175*u2)));
    B=u2/1024*(256+u2*(-128+u2*(74-47*u2)));
    ds=B*sins*(cos2sm+B/4*(coss*(-1+2*cos2sm^2)-B/6*cos2sm*(-3+4*sins^2)*(-3+4*cos2sm^2)));
    dist(i)=b*A*(s-ds)/1000;
    az(i)=mod(atan2(cosU2(i)*sin(lamda),cosU1(i)*sinU2(i)-sinU1(i)*cosU2(i)*cos(lamda))*R2D,360);
    baz(i)=mod(atan2(-cosU1(i)*sin(lamda),sinU1(i)*cosU2(i)-cosU1(i)*sinU2(i)*cos(lamda))*R2D,360);
    gcarc(i)=s*R2D;
end

end
