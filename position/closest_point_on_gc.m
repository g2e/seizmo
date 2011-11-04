function [lat3,lon3]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%CLOSEST_POINT_ON_GC    Return closest point on great circle to a point
%
%    Usage:    [lat3,lon3]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%    Description:
%     [LAT3,LON3]=CLOSEST_POINT_ON_GC(LAT1,LON1,LAT2,LON2,LAT,LON) finds
%     the point LAT3/LON3 on a great circle closest to the point(s) given
%     by LAT/LON.  The great circle(s) are given by LAT1/LON1 and
%     LAT2/LON2.  All LAT/LON must either be scalar or arrays with the same
%     number of elements.  This allows multiple points vs a great circle,
%     multiple great circles vs a point or 1 point per great circle.  All
%     inputs must be in degrees!
%
%    Notes:
%
%    Examples:
%     % Plot a great circle and its closest points to points of interest:
%     % random great circle
%     [lat1,lon1]=randlatlon;
%     [lat2,lon2]=randlatlon;
%     % random points
%     npts=round(10*rand);
%     [lat,lon]=randlatlon(npts);
%     % corresponding closest points
%     [lat3,lon3]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat,lon);
%     % draw world map
%     m_proj('hammer');
%     m_coast('patch',[0.6 1 0.6]);
%     % draw great circle
%     [gclat,gclon]=gc2latlon(lat1,lon1,lat2,lon2);
%     m_line(gclon',gclat','linewi',3,'color','b');
%     % draw gc arc from points to closest point on gc
%     [gclat,gclon]=gcarc2latlon(lat,lon,lat3,lon3);
%     gclon=unwrap(gclon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(gclon',gclat','linewi',3,'color','r');
%     % draw circles at points
%     m_line(lon3,lat3,'marker','o',...
%         'markerfacecolor','y','linestyle','none');
%     m_line(lon,lat,'marker','o',...
%         'markerfacecolor','m','linestyle','none');
%     % add grid to map
%     m_grid('xticklabels',[]);
%
%    See also: DEGDIST_FROM_GC, GC_INTERSECT, GC2LATLON, GCARC2LATLON,
%              GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Dec.  2, 2009 - improved example
%        June 26, 2010 - fix unwrap call in example
%        July 13, 2011 - code/doc cleanup
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 13, 2011 at 02:55 GMT

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
    error('seizmo:closest_point_on_gc:nonNumeric',...
        'All inputs must be numeric!');
end

% 1. convert A B C to xyz
A=[cosd(lon1).*cosd(lat1) sind(lon1).*cosd(lat1) sind(lat1)];
B=[cosd(lon2).*cosd(lat2) sind(lon2).*cosd(lat2) sind(lat2)];
C=[cosd(lon3).*cosd(lat3) sind(lon3).*cosd(lat3) sind(lat3)];

% 2. get perpendicular N to plane formed by gc AB
N=cross(A,B,2);

% 3. get perpendicular M to plane formed by gc NC
M=cross(N,C,2);

% 4. get perpendicular D to plane formed by gc NM
D=cross(M,N,2);

% 5. fix D if all points on gc AB are equidistant (choosing A)
%    - this may need to be tuned
bad=vecnorm(D,2)<10*eps;
if(any(bad)); D(bad,:)=A(bad,:); end

% 6. get D intersection with sphere
[lat1,lon1]=xyz2geocentric(D(:,1),D(:,2),D(:,3));
[lat2,lon2]=xyz2geocentric(-D(:,1),-D(:,2),-D(:,3));

% 7. get closer point along D
closer=sphericalinv(lat1,lon1,lat3,lon3)<sphericalinv(lat2,lon2,lat3,lon3);
if(any(closer))
    lat3(closer)=lat1(closer);
    lon3(closer)=lon1(closer);
end
closer=~closer;
if(any(closer))
    lat3(closer)=lat2(closer);
    lon3(closer)=lon2(closer);
end

% fix for single output
if(nargout<2)
    lat3=[lat3 lon3];
end

end
