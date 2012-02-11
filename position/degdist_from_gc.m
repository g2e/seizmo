function [dist]=degdist_from_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%DEGDIST_FROM_GC    Distance from a point on a sphere to a great circle
%
%    Usage:    dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%    Description:
%     DIST=DEGDIST_FROM_GC(LAT1,LON1,LAT2,LON2,LAT,LON) calculates the
%     shortest degree distance from a point on a sphere LAT/LON to a great
%     circle defined by the two points on a sphere LAT1/LON1 and LAT2/LON2.
%     All inputs must be in degrees!
%
%    Notes:
%     - Assumes positions are in geocentric coordinates.
%     - 180-DIST gives the farthest distance(s) to the great circle(s)
%
%    Examples:
%     % Find the distance from a series of
%     % random points to a random great circle:
%     [lat1,lon1]=randlatlon;
%     [lat2,lon2]=randlatlon;
%     npts=30;
%     [lat,lon]=randlatlon(npts);
%     dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon);
%     % Now plot the points and the great circle to visually confirm:
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [gclat,gclon]=gc2latlon(lat1,lon1,lat2,lon2);
%     gclon=unwrap(gclon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(gclon',gclat','linewi',3);
%     m_line(lon,lat,'linestyle','none','marker','o',...
%            'markerfacecolor','r','markersize',15);
%     m_text(lon,lat,num2str((1:npts)'),...
%            'hor','cen','ver','cap','color','w');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%     [(1:npts)' dist] % for confirmation
%
%    See also: CLOSEST_POINT_ON_GC, GC_INTERSECT, GC2LATLON, GCARC2LATLON,
%              GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%        July 13, 2011 - works now (forgot to normalize vectors)
%        Feb. 10, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

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
    error('seizmo:degdist_from_gc:nonNumeric',...
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

% 2. get unit vector N perpendicular to plane formed by gc AB
N=cross(A,B,2);
tmp=sum(abs(N).^2,2).^(1/2);
N=N./tmp(:,[1 1 1]);

% 3. get angle NOC abs(asind(NdotC))
dist=abs(asind(dot(C,N,2)));

% reshape to input
dist=reshape(dist,osz);

end
