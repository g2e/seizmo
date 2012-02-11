function [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2,npts)
%GC2LATLON    Returns points along great circle(s)
%
%    Usage:    [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2)
%              [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2,npts)
%
%    Description:
%     [LAT,LON]=GC2LATLON(LAT1,LON1,LAT2,LON2) calculates points LAT/LON
%     along the great circle(s) defined by LAT1/LON1 and LAT2/LON2 at
%     increments of 1 degree.  All LAT/LON arguments must either be scalar
%     or arrays with the same number of elements.  Outputs are returned so
%     that each row corresponds to a separate great circle given by the
%     input arguments.  All inputs must be in degrees!
%
%     [LAT,LON]=GC2LATLON(LAT1,LON1,LAT2,LON2,NPTS) specifies an
%     alternative number of points for the great circle discretization.
%     The default value for NPTS is 361 points between 0 & 360 degrees.
%
%    Notes:
%     - Assumes positions are given in geocentric coordinates.
%
%    Examples:
%     % Plot many great circles:
%     [lat1,lon1]=randlatlon(25);
%     [lat2,lon2]=randlatlon(25);
%     m_proj('robinson');
%     m_coast('color',[0 .6 0]);
%     [lat,lon]=gc2latlon(lat1,lon1,lat2,lon2);
%     lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
%     m_line(lon',lat','linewi',3);
%     m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: GCARC2LATLON, DEGDIST_FROM_GC, GC_INTERSECT,
%              GCARC_INTERSECT, CLOSEST_POINT_ON_GC

%     Version History:
%        Nov. 15, 2009 - initial version
%        Dec. 13, 2010 - fix wrap-around bug in example
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 10, 2012 - doc update, skip expansion
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% default/check npts
if(nargin==4 || isempty(npts)); npts=361; end
if(~isreal(npts) || ~isscalar(npts) || npts~=fix(npts) || npts<2)
    error('seizmo:gc2latlon:badSteps',...
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
    error('seizmo:gc2latlon:nonNumeric',...
        'All inputs must be real-valued arrays!');
elseif(sum(n~=1)>1 && numel(unique(n(n~=1)))>1)
    error('seizmo:gc2latlon:badSize',...
        'All inputs must have an equal number of elements or be scalar!');
end

% get great circles
[dist,az]=sphericalinv(lat1,lon1,lat2,lon2); n=numel(az);
dist=linspace(0,360,npts);
[lat,lon]=deal(nan(n,npts));
for i=1:n
    [lat(i,:),lon(i,:)]=sphericalfwd(lat1(i),lon1(i),dist,az(i));
end

end
