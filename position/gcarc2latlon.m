function [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2,npts)
%GCARC2LATLON    Returns points along great circle arc(s)
%
%    Usage:    [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2)
%              [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2,npts)
%
%    Description: [LAT,LON]=GCARC2LATLON(LAT1,LON1,LAT2,LON2) calculates
%     100 points LAT/LON along the great circle arc(s) defined by LAT1/LON1
%     and LAT2/LON2.  All LAT/LON arguments must either be scalar or arrays
%     with the same number of elements.  Outputs are returned so that each
%     row corresponds to a separate great circle given by the input
%     arguments.    All inputs must be in degrees!
%
%     [LAT,LON]=GCARC2LATLON(LAT1,LON1,LAT2,LON2,NPTS) specifies an
%     alternative number of points for the great circle arc discretization.
%     The default value for NPTS is 100 points.
%
%    Notes:
%
%    Examples:
%     Plot several great circle arcs:
%      lat1=90*2*(rand(5,1)-.5); lon1=180*2*(rand(5,1)-.5);
%      lat2=90*2*(rand(5,1)-.5); lon2=180*2*(rand(5,1)-.5);
%      m_proj('robinson');
%      m_coast('color',[0 .6 0]);
%      [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2);
%      m_line(lon',lat','linewi',3);
%      m_grid('linestyle','none','box','fancy','tickdir','out');
%
%    See also: GC2LATLON, DEGDIST_FROM_GC, GC_INTERSECT, GCARC_INTERSECT,
%              CLOSEST_POINT_ON_GC

%     Version History:
%        Nov. 15, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2009 at 19:35 GMT

% todo:

% check nargin
msg=nargchk(4,5,nargin);
if(~isempty(msg)); error(msg); end

% default/check npts
if(nargin==4 || isempty(npts)); npts=100; end
if(~isreal(npts) || ~isscalar(npts) || npts~=fix(npts) || npts<2)
    error('seizmo:gcarc2latlon:badSteps',...
        'STEPS must be a scalar integer >1!');
end

% lat/lon to column vector
lat1=lat1(:); lon1=lon1(:);
lat2=lat2(:); lon2=lon2(:);

% check lat/lon
sz1=size(lat1); sz2=size(lon1);
sz3=size(lat2); sz4=size(lon2);
n1=prod(sz1); n2=prod(sz2);
n3=prod(sz3); n4=prod(sz4);

% basic check inputs
if(~isreal(lat1) || ~isreal(lon1) || ~isreal(lat2) || ~isreal(lon2))
    error('seizmo:gcarc2latlon:nonNumeric',...
        'All inputs must be numeric!');
elseif(any([n1 n2 n3 n4]==0))
    lat=zeros(0,npts);
    lon=zeros(0,npts);
    return;
end

% expand scalars
if(n1==1); lat1=repmat(lat1,sz2); n1=n2; sz1=sz2; end
if(n2==1); lon1=repmat(lon1,sz1); sz2=sz1; end
if(n3==1); lat2=repmat(lat2,sz4); n3=n4; sz3=sz4; end
if(n4==1); lon2=repmat(lon2,sz3); sz4=sz3; end

% check inputs pairs
if(~isequal(sz1,sz2) || ~isequal(sz3,sz4))
    error('seizmo:gcarc2latlon:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% expand scalars
if(n1==1 && n3~=1); lat1=repmat(lat1,sz3); lon1=repmat(lon1,sz3); n1=n3; sz1=sz3; end
if(n3==1 && n1~=1); lat2=repmat(lat2,sz1); lon2=repmat(lon2,sz1); sz3=sz1; end

% cross check inputs
if(~isequal(sz1,sz3))
    error('seizmo:gcarc2latlon:nonscalarUnequalArrays',...
        'Input arrays need to be scalar or have equal size!');
end

% get great circle arcs
[dist,az]=sphericalinv(lat1,lon1,lat2,lon2);
lat=nan(n1,npts);
lon=nan(n1,npts);
for i=1:n1
    [lat(i,:),lon(i,:)]=sphericalfwd(lat1(i),lon1(i),...
        linspace(0,dist(i),npts),az(i));
end

end
