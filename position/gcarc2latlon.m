function [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2,npts)
%GCARC2LATLON    Returns points along great circle arc(s)
%
%    Usage:    [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2)
%              [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2,npts)
%
%    Description:
%     [LAT,LON]=GCARC2LATLON(LAT1,LON1,LAT2,LON2) calculates 100 points
%     LAT/LON along the great circle arc(s) defined by LAT1/LON1 and
%     LAT2/LON2.  All LAT/LON arguments must either be scalar or arrays
%     with the same number of elements.  Outputs are returned so that each
%     row corresponds to a separate great circle given by the input
%     arguments.    All inputs must be in degrees!
%
%     [LAT,LON]=GCARC2LATLON(LAT1,LON1,LAT2,LON2,NPTS) specifies an
%     alternative number of points for the great circle arc discretization.
%     The default value for NPTS is 100 points.
%
%    Notes:
%     - Assumes positions are given in geocentric coordinates.
%
%    Examples:
%     % Plot many great circle arcs:
%     [lat1,lon1]=randlatlon(100);
%     [lat2,lon2]=randlatlon(100);
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: GC2LATLON, DEGDIST_FROM_GC, GC_INTERSECT, GCARC_INTERSECT,
%              CLOSEST_POINT_ON_GC

%     Version History:
%        Nov. 15, 2009 - initial version
%        Dec. 13, 2010 - fix wrap-around bug in example
%        Jan. 22, 2011 - nargchk fix, code formatting
%        Feb. 10, 2012 - doc update, skip expansion
%        Oct. 15, 2012 - bugfix: expand only lat1/lon1 as needed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2012 at 19:35 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% default/check npts
if(nargin==4 || isempty(npts)); npts=100; end
if(~isreal(npts) || ~isscalar(npts) || npts~=fix(npts) || npts<2)
    error('seizmo:gcarc2latlon:badSteps',...
        'NPTS must be an integer >1!');
end

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);

% size up inputs
n(1)=numel(lat1); n(2)=numel(lon1);
n(3)=numel(lat2); n(4)=numel(lon2);

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2))
    error('seizmo:gcarc2latlon:nonNumeric',...
        'All inputs must be real-valued arrays!');
elseif(sum(n~=1)>1 && numel(unique(n(n~=1)))>1)
    error('seizmo:gcarc2latlon:badSize',...
        'All inputs must have an equal number of elements or be scalar!');
end

% get great circle arcs
[dist,az]=sphericalinv(lat1,lon1,lat2,lon2); n=numel(az);
if(isscalar(lat1)); lat1=lat1(ones(n,1)); end
if(isscalar(lon1)); lon1=lon1(ones(n,1)); end
[lat,lon]=deal(nan(n,npts));
for i=1:n
    [lat(i,:),lon(i,:)]=sphericalfwd(lat1(i),lon1(i),...
        linspace(0,dist(i),npts),az(i));
end

end
