function [lgc]=sph_poly_in(vlatlon,poly,plat,plon,overlap)
%SPH_POLY_IN    Returns indexing for which polygon points are in
%
%    Usage:    lgc=sph_poly_in(vlatlon,poly,plat,plon)
%              lgc=sph_poly_in(vlatlon,poly,plat,plon,overlap)
%
%    Description:
%     LGC=SPH_POLY_IN(VLATLON,POLY,PLAT,PLON) returns an indexing array LGC
%     indicating which polygon in POLY contains the points given by
%     PLAT/PLON.  Each point will only be found in one polygon (see the
%     next option to a point to be in multiple polygons).  VLATLON is the
%     vertex locations in latitude and longitude (N vertices as an Nx2
%     array of [LAT LON]) while POLY is the vertex indices in each polygon
%     (M polygons with Z vertices as an MxZ array).  PLAT/PLON must have
%     the same number of elements.  LGC is a NPTSxNPOLY logical array where
%     the each row is a logical vector indicating which
%     polygon the point is in.
%
%     LGC=SPH_POLY_IN(VLATLON,POLY,PLAT,PLON,OVERLAP) allows points in
%     multiple polygons.  The default (FALSE) gives only the first polygon
%     to contain each point.  Setting OVERLAP to TRUE allows points that
%     are in or on the edge of multiple polygons to return indexing
%     indicating all polygons that contain those points.
%
%    Notes:
%     - Polygons are always seen as the side with lesser area.
%
%    Examples:
%     % Which polygon is St. Louis in?
%     [v,tri]=sph_tri_init;
%     [lat,lon]=xyz2geocentric(v(:,1),v(:,2),v(:,3));
%     sph_poly_in([lat lon],tri,38.649,-90.305)
%
%    See also: SPH_POLY_AREA, SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_TRI_AUTO

%     Version History:
%        Nov. 17, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix, minor doc formatting
%        Feb. 18, 2011 - major revision: use vertex & index input/output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% check vlatlon
if(~isreal(vlatlon) || size(vlatlon,2)~=2)
    error('seizmo:sph_poly_in:badInput',...
        'VLATLON must be a real-valued Nx2 array!');
elseif(~isreal(poly) || size(poly,2)<3)
    error('seizmo:sph_poly_in:badInput',...
        'POLY must be a real-valued MxZ array where Z>=3!');
end

% check plat/plon (must be equal sized)
if(~isreal(plat) || ~isreal(plon))
    error('seizmo:sph_poly_in:badInput',...
        'PLAT & PLON must be real-valued arrays!');
end

% force column vector for point lat/lon
plat=plat(:); plon=plon(:);

% expand scalar point input
if(isscalar(plat)); plat=plat(ones(numel(plon),1)); end
if(isscalar(plon)); plon=plon(ones(numel(plat),1)); end
npts=numel(plat);

% require equal number of elements
if(npts~=numel(plon))
    error('seizmo:sph_poly_in:badInput',...
        'PLAT & PLON must be scalar or have the same number of elements!');
end

% default overlap
if(nargin==4 || isempty(overlap)); overlap=false; end

% check overlap
if(~islogical(overlap) || ~isscalar(overlap))
    error('seizmo:sph_poly_in:badInput',...
        'OVERLAP must be a scalar logical!');
end

% number of vertices/polygons
[npoly,nv]=size(poly);
nvtot=size(vlatlon,1);

% centers of polygons
[x,y,z]=geocentric2xyz(vlatlon(:,1),vlatlon(:,2),1);
[clat,clon]=xyz2geocentric(sum(x(poly),2),sum(y(poly),2),sum(z(poly),2));

% distances to vertices
dv=sphericalinv(clat(:,ones(1,nv)),clon(:,ones(1,nv)),...
    vlatlon(poly),vlatlon(nvtot+poly));

% distances to points from centers
dp=sphericalinv(clat(:,ones(1,npts))',clon(:,ones(1,npts))',...
    plat(:,ones(1,npoly)),plon(:,ones(1,npoly)));

% which polygons may contain the point
maxdv=max(dv,[],2);
maybe=dp<=maxdv(:,ones(1,npts))';

% get xyz values
[cx,cy,cz]=geocentric2xyz(clat,clon,1);
[px,py,pz]=geocentric2xyz(plat,plon,1);

% loop over each point
lgc=false(npts,npoly);
for i=1:npts
    for j=find(maybe(i,:))
        if(winding([cx(j) cy(j) cz(j)],...
                [x(poly(j,:)) y(poly(j,:)) z(poly(j,:))],...
                [px(i) py(i) pz(i)]))
            lgc(i,j)=true;
            if(~overlap); break; end
        end
    end
end

end

function [lgc]=winding(c,v,p)

% close polygon
v=[v; v(1,:)];

% project onto tangent plane
p=p/(c*p');
nv=c*v';
v=v./nv([1 1 1],:)';

% get vectors from point to vertices
pv=v-p(ones(1,size(v,1)),:);

% if p == v then is a vertex and we count it
if(any(~vecnorm(pv,2))); lgc=true; return; end

% get angle between vectors
normpv=vecnorm(pv,2);
costheta=sum(pv(1:end-1,:).*pv(2:end,:),2)...
    ./(normpv(1:end-1).*normpv(2:end));
ccw=sum(p(ones(1,size(pv,1)-1),:).*cross(pv(2:end,:),pv(1:end-1,:),2),2)<0;
theta=((-1).^ccw).*abs(acos(costheta));
w=sum(theta)/(2*pi);
lgc=round(w)~=0;

end
