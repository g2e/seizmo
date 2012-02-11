function [lat3,lon3]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%CLOSEST_POINT_ON_GC    Return closest point on great circle to a point
%
%    Usage:    [lat,lon]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%
%    Description:
%     [LAT,LON]=CLOSEST_POINT_ON_GC(LAT1,LON1,LAT2,LON2,LAT3,LON3) finds
%     the point(s) LAT/LON on great circle(s) closest to the point(s) in
%     LAT3/LON3.  The great circle(s) are given by the positions in
%     LAT1/LON1 and LAT2/LON2.  All inputs must be in degrees!
%
%    Notes:
%     - Assumes positions are in geocentric coordinates.
%     - The farthest points are those antipodal to the returned ones.
%
%    Examples:
%     % Plot a great circle and positions on it that are
%     % closest to some positions off the great circle:
%     [lat1,lon1]=randlatlon;
%     [lat2,lon2]=randlatlon;
%     [lat3,lon3]=randlatlon(100);
%     % corresponding closest points
%     [lat,lon]=closest_point_on_gc(lat1,lon1,lat2,lon2,lat3,lon3);
%     % draw world map
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     % draw great circle
%     [gclat,gclon]=gc2latlon(lat1,lon1,lat2,lon2);
%     gclon=unwrap(gclon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(gclon',gclat','linewi',3,'color','b');
%     % draw gc arc from points to closest points on gc
%     [gclat,gclon]=gcarc2latlon(lat3,lon3,lat,lon);
%     gclon=unwrap(gclon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(gclon',gclat','linewi',3,'color','r');
%     % draw circles at points
%     m_line(lon3,lat3,'marker','o',...
%         'markerfacecolor','y','linestyle','none');
%     m_line(lon,lat,'marker','o',...
%         'markerfacecolor','m','linestyle','none');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: DEGDIST_FROM_GC, GC_INTERSECT, GC2LATLON, GCARC2LATLON,
%              GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Dec.  2, 2009 - improved example
%        June 26, 2010 - fix unwrap call in example
%        July 13, 2011 - code/doc cleanup
%        Feb. 10, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 02:55 GMT

% todo:

% check nargin
error(nargchk(6,6,nargin));

% check lat/lon
sz{1}=size(lat1); sz{2}=size(lon1); sz{3}=size(lat2);
sz{4}=size(lon2); sz{5}=size(lat3); sz{6}=size(lon3);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3});
n(4)=prod(sz{4}); n(5)=prod(sz{5}); n(6)=prod(sz{6});

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) ...
        || ~isreal(lon2) || ~isreal(lat3) || ~isreal(lon3))
    error('seizmo:closest_point_on_gc:nonNumeric',...
        'All inputs must be real-valued arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:closest_point_on_gc:badSize',...
        'All inputs must be equal sized or scalar!');
end

% output size
if(sum(n~=1))
    osz=sz{find(n~=1,1)};
    on=prod(osz);
else
    osz=sz{1};
    on=1;
end

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);
lat3=lat3(:); lon3=lon3(:);

% expand scalars
if(n(1)==1); lat1=lat1(ones(on,1)); end
if(n(2)==1); lat1=lat1(ones(on,1)); end
if(n(3)==1); lat2=lat2(ones(on,1)); end
if(n(4)==1); lat2=lat2(ones(on,1)); end
if(n(5)==1); lat3=lat3(ones(on,1)); end
if(n(6)==1); lat3=lat3(ones(on,1)); end

% 1. get xyz
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
%    - this occurs if the point lies on vector perpendicular
%      to the plane formed by the great circle
%    - this may need to be tuned
bad=sqrt(sum(abs(D).^2,2))<10*sqrt(eps);
if(any(bad)); D(bad,:)=A(bad,:); end

% 6. get D intersection with sphere
[lat1,lon1]=xyz2geocentric(D(:,1),D(:,2),D(:,3));
[lat2,lon2]=xyz2geocentric(-D(:,1),-D(:,2),-D(:,3));

% 7. get closer point along D
two=haversine(lat1,lon1,lat3,lon3)>haversine(lat2,lon2,lat3,lon3);
[lat3,lon3]=deal(lat1,lon1);
[lat3(two),lon3(two)]=deal(lat2(two),lon2(two));

% reshape to input
lat3=reshape(lat3,osz); lon3=reshape(lon3,osz);

end
