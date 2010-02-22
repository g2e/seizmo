function [lat1,lon1,lat2,lon2]=gc_intersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%GC_INTERSECT    Return intersection points between great circles
%
%    Usage:     [ilat1,ilon1,ilat2,ilon2]=...
%                   gc_intersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%
%    Description: [ILAT1,ILON1,ILAT2,ILON2]=GC_INTERSECT(...
%     LAT1,LON1,LAT2,LON2,LAT3,LON3,LAT4,LON4) finds the intersection
%     points of great circles given by points LAT1/LON1 and LAT2/LON2 with
%     great circles given by points LAT3/LON3 and LAT4/LON4.  Great circles
%     either intersect twice or are equal.  ILAT1/ILON1 gives one
%     intersection point and ILAT2/ILON2 gives the other (antipodal to the
%     first).  When two great circles are equal both intersection points
%     are set to NaN.  All LAT/LON must either be scalar or arrays with the
%     same number of elements.  This allows finding intersections between
%     one great circle and several others or to find intersections between
%     distinct pairs.  All inputs must be in degrees!
%
%    Notes:
%
%    Examples:
%     Plot 2 great circles and their intersections:
%      lat1=90*2*(rand(2,1)-.5); lon1=180*2*(rand(2,1)-.5);
%      lat2=90*2*(rand(2,1)-.5); lon2=180*2*(rand(2,1)-.5);
%      m_proj('robinson');
%      m_coast('color',[0 .6 0]);
%      [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2);
%      m_line(lon',lat','linewi',3);
%      [ilat1,ilon1,ilat2,ilon2]=gc_intersect(lat1(1),lon1(1),...
%          lat2(1),lon2(1),lat1(2),lon1(2),lat2(2),lon2(2));
%      m_line(ilon1,ilat1,'marker','o','markerfacecolor','y');
%      m_line(ilon2,ilat2,'marker','o','markerfacecolor','r');
%      m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
%              GCARC2LATLON, GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2009 at 19:35 GMT

% todo:

% check nargin
msg=nargchk(8,8,nargin);
if(~isempty(msg)); error(msg); end

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);
lat3=lat3(:); lon3=lon3(:);
lat4=lat4(:); lon4=lon4(:);

% check lat/lon
sz1=size(lat1); sz2=size(lon1);
sz3=size(lat2); sz4=size(lon2);
sz5=size(lat3); sz6=size(lon3);
sz7=size(lat4); sz8=size(lon4);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);
n5=prod(sz5); n6=prod(sz6);
n7=prod(sz7); n8=prod(sz8);

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2) || ...
        ~isreal(lat3) || ~isreal(lon3) || ~isreal(lat4) || ~isreal(lon4))
    error('seizmo:gc_intersect:nonNumeric',...
        'All inputs must be numeric!');
elseif(any([n1 n2 n3 n4 n5 n6 n7 n8]==0))
    lat1=zeros(0,1);
    lon1=zeros(0,1);
    lat2=zeros(0,1);
    lon2=zeros(0,1);
    return;
end

% expand scalars
if(n1==1); lat1=repmat(lat1,sz2); n1=n2; sz1=sz2; end
if(n2==1); lon1=repmat(lon1,sz1); sz2=sz1; end
if(n3==1); lat2=repmat(lat2,sz4); n3=n4; sz3=sz4; end
if(n4==1); lon2=repmat(lon2,sz3); sz4=sz3; end
if(n5==1); lat3=repmat(lat3,sz6); n5=n6; sz5=sz6; end
if(n6==1); lon3=repmat(lon3,sz5); sz6=sz5; end
if(n7==1); lat4=repmat(lat4,sz8); n7=n8; sz7=sz8; end
if(n8==1); lon4=repmat(lon4,sz7); sz8=sz7; end

% check input pairs
if(~isequal(sz1,sz2) || ~isequal(sz3,sz4) || ...
        ~isequal(sz5,sz6) || ~isequal(sz7,sz8))
    error('seizmo:gc_intersect:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n1==1 && n3~=1); lat1=repmat(lat1,sz3); lon1=repmat(lon1,sz3); n1=n3; sz1=sz3; end
if(n1==1 && n5~=1); lat1=repmat(lat1,sz5); lon1=repmat(lon1,sz5); n1=n5; sz1=sz5; end
if(n1==1 && n7~=1); lat1=repmat(lat1,sz7); lon1=repmat(lon1,sz7); n1=n7; sz1=sz7; end
if(n3==1 && n1~=1); lat2=repmat(lat2,sz1); lon2=repmat(lon2,sz1); n3=n1; sz3=sz1; end
if(n3==1 && n5~=1); lat2=repmat(lat2,sz5); lon2=repmat(lon2,sz5); n3=n5; sz3=sz5; end
if(n3==1 && n7~=1); lat2=repmat(lat2,sz7); lon2=repmat(lon2,sz7); n3=n7; sz3=sz7; end
if(n5==1 && n1~=1); lat3=repmat(lat3,sz1); lon3=repmat(lon3,sz1); n5=n1; sz5=sz1; end
if(n5==1 && n3~=1); lat3=repmat(lat3,sz3); lon3=repmat(lon3,sz3); n5=n3; sz5=sz3; end
if(n5==1 && n7~=1); lat3=repmat(lat3,sz7); lon3=repmat(lon3,sz7); n5=n7; sz5=sz7; end
if(n7==1 && n1~=1); lat4=repmat(lat4,sz1); lon4=repmat(lon4,sz1); n7=n1; sz7=sz1; end
if(n7==1 && n3~=1); lat4=repmat(lat4,sz3); lon4=repmat(lon4,sz3); n7=n3; sz7=sz3; end
if(n7==1 && n5~=1); lat4=repmat(lat4,sz5); lon4=repmat(lon4,sz5); sz7=sz5; end

% cross check inputs
if(~isequal(sz1,sz3,sz5,sz7))
    error('seizmo:gc_intersect:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% 1. convert A B C to xyz
A=[cosd(lon1).*cosd(lat1) sind(lon1).*cosd(lat1) sind(lat1)];
B=[cosd(lon2).*cosd(lat2) sind(lon2).*cosd(lat2) sind(lat2)];
C=[cosd(lon3).*cosd(lat3) sind(lon3).*cosd(lat3) sind(lat3)];
D=[cosd(lon4).*cosd(lat4) sind(lon4).*cosd(lat4) sind(lat4)];

% 2. get perpendicular M to plane formed by gc AB
M=cross(A,B,2);

% 3. get perpendicular N to plane formed by gc CD
N=cross(C,D,2);

% 4. get perpendicular E to plane formed by gc MN
E=cross(M,N,2);

% 5. get E intersection with sphere
[lat1,lon1]=xyz2geocentric(E(:,1),E(:,2),E(:,3));
[lat2,lon2]=xyz2geocentric(-E(:,1),-E(:,2),-E(:,3));

% 6. nans if E is basically zero - this may need to be tuned
bad=vecnorm(E,2)<10*eps;
if(any(bad))
    lat1(bad)=nan; lon1(bad)=nan;
    lat2(bad)=nan; lon2(bad)=nan;
    return;
end

end
