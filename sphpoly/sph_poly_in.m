function [lgc]=sph_poly_in(vlat,vlon,plat,plon,overlap)
%SPH_POLY_IN    Returns indexing for which polygon points are in
%
%    Usage:    lgc=sph_poly_in(vlat,vlon,plat,plon)
%              lgc=sph_poly_in(vlat,vlon,plat,plon,overlap)
%
%    Description:
%     LGC=SPH_POLY_IN(VLAT,VLON,PLAT,PLON) returns an indexing array LGC
%     indicating which polygon in VLAT/VLON contains the points given by
%     PLAT/PLON.  Each point will only be found in one polygon (see the
%     next option to a point to be in multiple polygons).  VLAT/VLON must
%     be equal sized and each column gives vertices for separate polygons.
%     PLAT/PLON must have the same number of elements.  LGC is a NPTSxNPOLY
%     logical array where the each row is a logical vector indicating which
%     polygon the point is in.
%
%     LGC=SPH_POLY_IN(VLAT,VLON,PLAT,PLON,OVERLAP) allows points to be in
%     multiple polygons.  The default (FALSE) gives only the first polygon
%     to contain each point.  Setting OVERLAP to TRUE allows points that
%     are in or on the edge of multiple polygons to return indexing
%     indicating all polygons that contain those points.
%
%    Notes:
%
%    Examples:
%     % Which polygon is St. Louis in?
%     [a,b,c]=sph_tri_init;
%     [lat,lon]=sph_tri_latlon(a,b,c);
%     sph_poly_in(lat,lon,38.649,-90.305)
%
%    See also: SPH_POLY_AREA, SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_TRI_AUTO,
%              SPH_TRI_LATLON, SPH_TRI_VERTICES

%     Version History:
%        Nov. 17, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix, minor doc formatting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% check vlat/vlon (must be equal sized)
if(~isreal(vlat) || ~isreal(vlon) || ~isequal(size(vlat),size(vlon)))
    error('seizmo:sph_poly_in:badInput',...
        'VLAT & VLON must be equal sized!');
end

% convert row vector to column
if(isvector(vlat)); vlat=vlat(:); end
if(isvector(vlon)); vlat=vlon(:); end

% force column vector for point lat/lon
plat=plat(:); plon=plon(:);
npts=size(plat,1);

% check plat/plon (must be equal sized)
if(~isreal(plat) || ~isreal(plon) || ~isequal(size(plat),size(plon)))
    error('seizmo:sph_poly_in:badInput',...
        'PLAT & PLON must be equal sized!');
end

% default overlap
if(nargin==4 || isempty(overlap)); overlap=false; end

% check overlap
if(~islogical(overlap) || ~isscalar(overlap))
    error('seizmo:sph_poly_in:badInput',...
        'OVERLAP must be a scalar logical!');
end

% number of vertices
[nv,npoly]=size(vlat);

% centers of polygons
[x,y,z]=geocentric2xyz(vlat,vlon,1);
[clat,clon]=xyz2geocentric(sum(x,1),sum(y,1),sum(z,1));

% distances to vertices
dv=sphericalinv(...
    submat(clat,1,ones(nv,1)),submat(clon,1,ones(nv,1)),vlat,vlon);

% distances to points
dp=sphericalinv(...
    submat(clat,1,ones(npts,1)),submat(clon,1,ones(npts,1)),...
    plat(:,ones(1,npoly)),plon(:,ones(1,npoly)));

% which polygons may contain the point
maxdv=max(dv);
ok=dp<=submat(maxdv,1,ones(npts,1));

% get xyz values
[cx,cy,cz]=geocentric2xyz(clat,clon,1);
[vx,vy,vz]=geocentric2xyz(vlat,vlon,1);
[px,py,pz]=geocentric2xyz(plat,plon,1);

% loop over each point
lgc=false(npts,npoly);
for i=1:npts
    oki=ok(i,:);
    for j=find(oki)
        if(winding([cx(j) cy(j) cz(j)],...
                [vx(:,j) vy(:,j) vz(:,j)],...
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
