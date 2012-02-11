function [lat,lon]=gcarc_intersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%GCARC_INTERSECT    Return intersection points between great circle arcs
%
%    Usage:     [lat,lon]=gcarc_intersect(...
%                   lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%
%    Description:
%     [LAT,LON]=GCARC_INTERSECT(LAT1,LON1,LAT2,LON2,LAT3,LON3,LAT4,LON4)
%     finds the intersection points of great circle arcs given by points
%     LAT1/LON1 and LAT2/LON2 with great circle arcs given by points
%     LAT3/LON3 and LAT4/LON4.  Arcs are the shortest path between points.
%     Great circle arcs either intersect once or not at all (arcs along the
%     same great circle are always treated as non-intersecting).  When two
%     great circle arcs do not intersect the LAT/LON for that pair is set
%     to NaN.  All inputs must be in degrees!
%
%    Notes:
%     - Assumes positions are in geocentric coordinates.
%
%    Examples:
%     % Plot several great circle arcs and their intersections:
%     [lat1,lon1]=randlatlon(25); [lat2,lon2]=randlatlon(25);
%     [lat3,lon3]=randlatlon(25); [lat4,lon4]=randlatlon(25);
%     % plot arcs
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     [lat,lon]=gcarc2latlon(lat3,lon3,lat4,lon4);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     % get intersections
%     [ilat1,ilon1]=gcarc_intersect(lat1,lon1,...
%         lat2,lon2,lat3,lon3,lat4,lon4);
%     % plot intersections
%     m_line(ilon1,ilat1,'marker','o',...
%         'markerfacecolor','y','linestyle','none');
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
%              GCARC2LATLON, GC_INTERSECT

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

% get intersection points of great circles
[flat,flon,slat,slon]=gc_intersect(...
    lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4);

% get great circle arc distances
dist1=haversine(lat1,lon1,lat2,lon2);
dist2=haversine(lat3,lon3,lat4,lon4);

% get distances of intersection points to arc end points
fdist1=haversine(lat1,lon1,flat,flon);
fdist2=haversine(lat2,lon2,flat,flon);
fdist3=haversine(lat3,lon3,flat,flon);
fdist4=haversine(lat4,lon4,flat,flon);
sdist1=haversine(lat1,lon1,slat,slon);
sdist2=haversine(lat2,lon2,slat,slon);
sdist3=haversine(lat3,lon3,slat,slon);
sdist4=haversine(lat4,lon4,slat,slon);

% allocate output
lat=nan(size(flat));
lon=nan(size(flon));

% check if distance from points to intersection adds up
% to the total arc length (or really, really close) as
% a method to determine if the intersection is on the arc
fon12=abs(fdist1+fdist2-dist1)<(10*dist1*sqrt(eps));
fon34=abs(fdist3+fdist4-dist2)<(10*dist2*sqrt(eps));
son12=abs(sdist1+sdist2-dist1)<(10*dist1*sqrt(eps));
son34=abs(sdist3+sdist4-dist2)<(10*dist2*sqrt(eps));

% intersection of circles needs to be
% on both arcs to be an arc intersection
fok=fon12 & fon34;
sok=son12 & son34;
if(any(fok))
    lat(fok)=flat(fok);
    lon(fok)=flon(fok);
end
if(any(sok))
    lat(sok)=slat(sok);
    lon(sok)=slon(sok);
end

end
