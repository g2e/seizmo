function [dist]=degdist_from_gc(lat1,lon1,lat2,lon2,lat3,lon3)
%DEGDIST_FROM_GC    Distance from a point on a sphere to a great circle
%
%    Usage:    dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%    Description: DIST=DEGDIST_FROM_GC(LAT1,LON1,LAT2,LON2,LAT,LON)
%     calculates the smallest degree distance from a point on a sphere
%     LAT/LON to a great circle defined by the two points on a sphere
%     LAT1/LON1 and LAT2/LON2.  All LAT/LON must either be scalar or arrays
%     with the same number of elements.  This allows operations like
%     finding the distance from multiple points to a single great circle or
%     finding the distance from a single point to several great circles.
%     DIST is always returned as a column vector.  All inputs must be in
%     degrees!
%
%    Notes:
%
%    Examples:
%     Find the distance from a series of
%     random points to a random great circle:
%      lat1=90*2*(rand-.5);
%      lon1=180*2*(rand-.5);
%      lat2=90*2*(rand-.5);
%      lon2=180*2*(rand-.5);
%      npts=round(10*rand);
%      lat=90*2*(rand(npts,1)-.5);
%      lon=180*2*(rand(npts,1)-.5);
%      dist=degdist_from_gc(lat1,lon1,lat2,lon2,lat,lon)
%
%     Now plot the points and the great circle to visually confirm:
%      m_proj('robinson');
%      m_coast('color',[0 .6 0]);
%      [gclat,gclon]=gc2latlon(lat1,lon1,lat2,lon2);
%      m_line(gclon',gclat','linewi',3);
%      m_line(lon,lat,'linestyle','none','marker','o',...
%             'markerfacecolor','r','markersize',10);
%      m_text(lon,lat,num2str((1:npts)'),'hor','cen','ver','cap');
%      m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: CLOSEST_POINT_ON_GC, GC_INTERSECT, GC2LATLON, GCARC2LATLON,
%              GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2009 at 07:15 GMT

% todo:

% check nargin
msg=nargchk(6,6,nargin);
if(~isempty(msg)); error(msg); end

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);
lat3=lat3(:); lon3=lon3(:);

% check lat/lon
sz1=size(lat1); sz2=size(lon1);
sz3=size(lat2); sz4=size(lon2);
sz5=size(lat3); sz6=size(lon3);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);
n5=prod(sz5); n6=prod(sz6);

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2) || ...
        ~isreal(lat3) || ~isreal(lon3))
    error('seizmo:degdist_from_gc:nonNumeric',...
        'All inputs must be numeric!');
elseif(any([n1 n2 n3 n4 n5 n6]==0))
    dist=zeros(0,1);
    return;
end

% expand scalars
if(n1==1); lat1=repmat(lat1,sz2); n1=n2; sz1=sz2; end
if(n2==1); lon1=repmat(lon1,sz1); sz2=sz1; end
if(n3==1); lat2=repmat(lat2,sz4); n3=n4; sz3=sz4; end
if(n4==1); lon2=repmat(lon2,sz3); sz4=sz3; end
if(n5==1); lat3=repmat(lat3,sz6); n5=n6; sz5=sz6; end
if(n6==1); lon3=repmat(lon3,sz5); sz6=sz5; end

% check inputs pairs
if(~isequal(sz1,sz2) || ~isequal(sz3,sz4) || ~isequal(sz5,sz6))
    error('seizmo:degdist_from_gc:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n1==1 && n3~=1); lat1=repmat(lat1,sz3); lon1=repmat(lon1,sz3); n1=n3; sz1=sz3; end
if(n1==1 && n5~=1); lat1=repmat(lat1,sz5); lon1=repmat(lon1,sz5); n1=n5; sz1=sz5; end
if(n3==1 && n1~=1); lat2=repmat(lat2,sz1); lon2=repmat(lon2,sz1); n3=n1; sz3=sz1; end
if(n3==1 && n5~=1); lat2=repmat(lat2,sz5); lon2=repmat(lon2,sz5); n3=n5; sz3=sz5; end
if(n5==1 && n1~=1); lat3=repmat(lat3,sz1); lon3=repmat(lon3,sz1); n5=n1; sz5=sz1; end
if(n5==1 && n3~=1); lat3=repmat(lat3,sz3); lon3=repmat(lon3,sz3); sz5=sz3; end

% cross check inputs
if(~isequal(sz1,sz3,sz5))
    error('seizmo:degdist_from_gc:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% 1. convert A B C to xyz
A=[cosd(lon1).*cosd(lat1) sind(lon1).*cosd(lat1) sind(lat1)];
B=[cosd(lon2).*cosd(lat2) sind(lon2).*cosd(lat2) sind(lat2)];
C=[cosd(lon3).*cosd(lat3) sind(lon3).*cosd(lat3) sind(lat3)];

% 2. get unit vector N perpendicular to plane formed by gc AB
N=cross(A,B,2);

% 3. get angle NOC abs(asind(NdotC))
dist=abs(asind(dot(C,N,2)));

end
