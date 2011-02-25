function [v,tri]=sph_tri_auto(maxdeg)
%SPH_TRI_AUTO    Automatically create spherical triangle polyhedra
%
%    Usage:    [v,tri]=sph_tri_auto(maxdeg)
%
%    Description:
%     [V,TRI]=SPH_TRI_AUTO(MAXDEG) automatically creates a spherical
%     triangle polyhedra where the vertices are all approximately less than
%     MAXDEG degrees apart.  Distortion due to projection of the triangle
%     vertices onto the unit sphere may lead to some vertices exceeding
%     MAXDEG.
%
%    Notes:
%     - Starts with an icosahedron via SPH_TRI_INIT
%     - The nearest vertices on an icosahedron are 63.435 deg
%
%    Examples:
%     % Make a polyhedra for Earth requiring ~50km spacing between vertices
%     % (warning this makes a >400000 faceted polygon!):
%     [v,tri]=sph_tri_auto(50*180/(6371*pi));
%     v=v*6371;
%
%     % Polyhedra with ~4deg vertex spacing:
%     [v,tri]=sph_tri_auto(4);
%     tri=tri(:,[1:3 1]); % close triangles
%     nv=size(v,1); % number of vertices
%     patch(v(tri)',v(nv+tri)',v(2*nv+tri)','r');
%     axis vis3d
%     axis off
%     rotate3d on
%
%    See also: SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov. 12, 2009 - initial version
%        Nov. 17, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix, minor doc formmatting
%        Feb. 18, 2011 - major revision: use vertex & index output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check maxdeg
if(~isscalar(maxdeg) || ~isreal(maxdeg) || maxdeg<=0)
    error('seizmo:sph_tri_auto:badInput',...
        'MAXDEG must be a positive real-valued scalar (in degrees)!');
end

% vertex spacing
% faces vert degdist
% 4     4    109.4712206344907
% 8     6    90
% 20    12   63.434948822922010
vd=63.434948822922010;

% splitting order
s=ceil(vd/maxdeg);

% factored splitting order
% - less distortion if we take small steps
s=factor(s);

% spherical triangle polyhedra
[v,tri]=sph_tri_init;
[v,tri]=sph_tri_split(v,tri,numel(s),s);

end
