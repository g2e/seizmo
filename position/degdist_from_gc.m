function [dist]=degdist_from_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%DEGDIST_FROM_GC    Distance from a point on a sphere to a great circle
%
%    Usage:    dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%    Description:
%     DIST=DEGDIST_FROM_GC(LAT1,LON1,LAT2,LON2,LAT,LON) calculates the
%     smallest degree distance from a point on a sphere LAT/LON to a great
%     circle defined by the two points on a sphere LAT1/LON1 and LAT2/LON2.
%     All LAT/LON must either be scalar or arrays with the same number of
%     elements.  This allows operations like finding the distance from
%     multiple points to a single great circle or finding the distance from
%     a single point to several great circles.  DIST is always returned as
%     a column vector.  All inputs must be in degrees!
%
%    Notes:
%
%    Examples:
%     % Find the distance from a series of
%     % random points to a random great circle:
%     [lat1,lon1]=randlatlon;
%     [lat2,lon2]=randlatlon;
%     npts=round(10*rand);
%     [lat,lon]=randlatlon(npts);
%     dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%     % Now plot the points and the great circle to visually confirm:
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [gclat,gclon]=gc2latlon(lat1,lon1,lat2,lon2);
%     m_line(gclon',gclat','linewi',3);
%     m_line(lon,lat,'linestyle','none','marker','o',...
%            'markerfacecolor','r','markersize',10);
%     m_text(lon,lat,num2str((1:npts)'),'hor','cen','ver','cap');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: CLOSEST_POINT_ON_GC, GC_INTERSECT, GC2LATLON, GCARC2LATLON,
%              GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%        July 13, 2011 - works now (forgot to normalize vectors)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 13, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(6,6,nargin));

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);
lat3=lat3(:); lon3=lon3(:);

% expand scalars
[lat1,lon1,lat2,lon2,lat3,lon3]=expandscalars(...
    lat1,lon1,lat2,lon2,lat3,lon3);

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2) || ...
        ~isreal(lat3) || ~isreal(lon3))
    error('seizmo:degdist_from_gc:nonNumeric',...
        'All inputs must be numeric!');
end

% 1. convert A B C to xyz (assumes a sphere)
A=[cosd(lon1).*cosd(lat1) sind(lon1).*cosd(lat1) sind(lat1)];
B=[cosd(lon2).*cosd(lat2) sind(lon2).*cosd(lat2) sind(lat2)];
C=[cosd(lon3).*cosd(lat3) sind(lon3).*cosd(lat3) sind(lat3)];

% 2. get unit vector N perpendicular to plane formed by gc AB
N=cross(A,B,2);
tmp=sum(abs(N).^2,2).^(1/2);
N=N./tmp(:,[1 1 1]);

% 3. get angle NOC abs(asind(NdotC))
dist=abs(asind(dot(C,N,2)));

end
