function [lat,lon]=gcarc_intersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%GCARC_INTERSECT    Return intersection points between great circle arcs
%
%    Usage:     [lat,lon]=gcarc_intersect(...
%                   lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%
%    Description: [LAT,LON]=GCARC_INTERSECT(LAT1,LON1,LAT2,LON2,...
%     LAT3,LON3,LAT4,LON4) finds the intersection points of great circle
%     arcs given by points LAT1/LON1 and LAT2/LON2 with great circle arcs
%     given by points LAT3/LON3 and LAT4/LON4.  Arcs are the short arc
%     between points.  Great circle arcs either intersect once or not at
%     all (arcs along the same great circle are always treated as non-
%     intersecting).  When two great circle arcs do not intersect the
%     return LAT/LON for that pair is set to NaN.  All LAT/LON inputs must
%     either be scalar or arrays with the same number of elements.  This
%     allows finding intersections between one great circle arc and several
%     others or to find intersections between distinct pairs.  All inputs
%     must be in degrees!
%
%    Notes:
%
%    Examples:
%     Plot several great circle arcs and their intersections:
%      lat1=90*2*(rand(5,1)-.5); lon1=180*2*(rand(5,1)-.5);
%      lat2=90*2*(rand(5,1)-.5); lon2=180*2*(rand(5,1)-.5);
%      lat3=90*2*(rand(5,1)-.5); lon3=180*2*(rand(5,1)-.5);
%      lat4=90*2*(rand(5,1)-.5); lon4=180*2*(rand(5,1)-.5);
%      m_proj('robinson');
%      m_coast('color',[0 .6 0]);
%      [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2);
%      m_line(lon',lat','linewi',3);
%      [lat,lon]=gcarc2latlon(lat3,lon3,lat4,lon4);
%      m_line(lon',lat','linewi',3);
%      [ilat1,ilon1]=gcarc_intersect(lat1,lon1,...
%          lat2,lon2,lat3,lon3,lat4,lon4);
%      m_line(ilon1,ilat1,'marker','o','markerfacecolor','y');
%      m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
%              GCARC2LATLON, GC_INTERSECT

%     Version History:
%        Nov. 15, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(8,8,nargin));

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);
lat3=lat3(:); lon3=lon3(:);
lat4=lat4(:); lon4=lon4(:);

% get intersection points of great circles
[flat,flon,slat,slon]=...
    gc_intersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4);

% get great circle arc distances
dist1=sphericalinv(lat1,lon1,lat2,lon2);
dist2=sphericalinv(lat3,lon3,lat4,lon4);

% get distances of intersection points to arc end points
fdist1=sphericalinv(lat1,lon1,flat,flon);
fdist2=sphericalinv(lat2,lon2,flat,flon);
fdist3=sphericalinv(lat3,lon3,flat,flon);
fdist4=sphericalinv(lat4,lon4,flat,flon);
sdist1=sphericalinv(lat1,lon1,slat,slon);
sdist2=sphericalinv(lat2,lon2,slat,slon);
sdist3=sphericalinv(lat3,lon3,slat,slon);
sdist4=sphericalinv(lat4,lon4,slat,slon);

% allocate output
lat=nan(size(flat));
lon=nan(size(flon));

% check if between - this may need to be tuned
fon12=(fdist1+fdist2-dist1)<(10*sqrt(eps)*dist1);
fon34=(fdist3+fdist4-dist2)<(10*sqrt(eps)*dist2);
son12=(sdist1+sdist2-dist1)<(10*sqrt(eps)*dist1);
son34=(sdist3+sdist4-dist2)<(10*sqrt(eps)*dist2);

% assign acceptable to output
fok=fon12&fon34;
sok=son12&son34;
if(any(fok))
    lat(fok)=flat(fok);
    lon(fok)=flon(fok);
end
if(any(sok))
    lat(sok)=slat(sok);
    lon(sok)=slon(sok);
end

end
