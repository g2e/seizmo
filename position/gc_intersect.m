function [lat1,lon1,lat2,lon2]=gc_intersect(...
    lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%GC_INTERSECT    Return intersection points between great circles
%
%    Usage:     [ilat1,ilon1,ilat2,ilon2]=gc_intersect(...
%                lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%
%    Description:
%     [ILAT1,ILON1,ILAT2,ILON2]=GC_INTERSECT(...
%         LAT1,LON1,LAT2,LON2,LAT3,LON3,LAT4,LON4) finds the intersection
%     points of great circles given by points LAT1/LON1 and LAT2/LON2 with
%     great circles given by points LAT3/LON3 and LAT4/LON4.  Great circles
%     either intersect twice or are equal.  ILAT1/ILON1 has 1 intersection
%     point and ILAT2/ILON2 gives the other (antipodal to the first).  When
%     two great circles are equal both intersection points are set to NaNs.
%     All LAT/LON must either be scalar or arrays with the same number of
%     elements.  This allows finding intersections between one great circle
%     and several others or to find intersections between distinct pairs.
%     All inputs must be in degrees!
%
%    Notes:
%     - Assumes positions are in geocentric coordinates.
%
%    Examples:
%     % Plot great circles and their intersections:
%     [lat1,lon1]=randlatlon(1); [lat2,lon2]=randlatlon(1);
%     [lat3,lon3]=randlatlon(5); [lat4,lon4]=randlatlon(5);
%     % plot circles
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     [lat,lon]=gc2latlon(lat3,lon3,lat4,lon4);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     % get intersections
%     [ilat1,ilon1,ilat2,ilon2]=gc_intersect(lat1,lon1,...
%         lat2,lon2,lat3,lon3,lat4,lon4);
%     % plot intersections
%     m_line(ilon1,ilat1,'marker','o',...
%         'markerfacecolor','y','linestyle','none');
%     m_line(ilon2,ilat2,'marker','o',...
%         'markerfacecolor','r','linestyle','none');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
%              GCARC2LATLON, GCARC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(8,8,nargin));

% check lat/lon
sz{1}=size(lat1); sz{2}=size(lon1); sz{3}=size(lat2); sz{4}=size(lon2);
sz{5}=size(lat3); sz{6}=size(lon3); sz{7}=size(lat4); sz{8}=size(lon4);
n(1)=prod(sz{1}); n(2)=prod(sz{2}); n(3)=prod(sz{3}); n(4)=prod(sz{4});
n(5)=prod(sz{5}); n(6)=prod(sz{6}); n(7)=prod(sz{7}); n(8)=prod(sz{8});

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2) || ...
        ~isreal(lat3) || ~isreal(lon3) || ~isreal(lat4) || ~isreal(lon4))
    error('seizmo:gc_intersect:nonNumeric',...
        'All inputs must be real-valued arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:gc_intersect:badSize',...
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
lat4=lat4(:); lon4=lon4(:);

% expand scalars
if(n(1)==1); lat1=lat1(ones(on,1)); end
if(n(2)==1); lat1=lat1(ones(on,1)); end
if(n(3)==1); lat2=lat2(ones(on,1)); end
if(n(4)==1); lat2=lat2(ones(on,1)); end
if(n(5)==1); lat3=lat3(ones(on,1)); end
if(n(6)==1); lat3=lat3(ones(on,1)); end
if(n(7)==1); lat4=lat4(ones(on,1)); end
if(n(8)==1); lat4=lat4(ones(on,1)); end

% 1. get xyz
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

% 5. get E intersections with sphere
[lat1,lon1]=xyz2geocentric(E(:,1),E(:,2),E(:,3));
[lat2,lon2]=xyz2geocentric(-E(:,1),-E(:,2),-E(:,3));

% 6. nans if E is basically zero - this may need to be tuned
bad=sqrt(sum(abs(E).^2,2))<10*sqrt(eps);
if(any(bad))
    lat1(bad)=nan; lon1(bad)=nan;
    lat2(bad)=nan; lon2(bad)=nan;
    return;
end

% reshape to input
lat1=reshape(lat1,osz); lon1=reshape(lon1,osz);
lat2=reshape(lat2,osz); lon2=reshape(lon2,osz);

end
