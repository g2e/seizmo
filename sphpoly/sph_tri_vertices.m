function [vertices]=sph_tri_vertices(a,b,c)
%SPH_TRI_VERTICES    Unique vertices of a spherical triangle polyhedra
%
%    Usage:    vertices=sph_tri_vertices(a,b,c)
%
%    Description:
%     VERTICES=SPH_TRI_VERTICES(A,B,C) extracts the unique vertices from
%     the spherical triangle polyhedra formed by A, B, & C.  VERTICES is an
%     Mx3 array of xyz values where M is the number of vertices.  If A, B,
%     & C are Nx3 arrays (where N is the number of faces of the polyhedra)
%     then M should be N/2+2.
%
%    Notes:
%     - If you get a warning, this means that a single vertex is likely
%       being seen as multiple vertices.  You will need to remove any
%       vertices that are essentially at the same location (leaving one).
%
%    Examples:
%     % Get a grid of points in latitude and longitude that are all
%     % approximately equally spaced at roughly 1deg:
%     [a,b,c]=sph_tri_auto(1);
%     v=sph_tri_vertices(a,b,c);
%     [lat,lon]=xyz2geocentric(v(:,1),v(:,2),v(:,3));
%
%    See also: SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_TRI_AUTO, SPH_TRI_LATLON,
%              SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov. 12, 2009 - initial version
%        Nov. 17, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix, warn fix, minor doc formmatting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% check arguments
sa=size(a); sb=size(b); sc=size(c);
if(~isequal(sa,sb,sc) || numel([sa sb sc])>6 ...
        || sa(2)~=3 || sb(2)~=3 || sc(2)~=3)
    error('seizmo:sph_tri_vertices:badInput',...
        'A, B, & C must be equal-sized Mx3 arrays!');
elseif(any(abs(vecnorm(cat(1,a,b,c),2)-1)>sqrt(3*eps*eps)))
    error('seizmo:sph_tri_vertices:badInput',...
        'A, B, & C vertices must be on the unit sphere!');
end

% let unique handle sorting and elimination
vertices=unique(cat(1,a,b,c),'rows');

% now lets check
if(size(vertices,1)~=(sa/2+2))
    % something is wrong
    warning('seizmo:sph_tri_vertices:verticesCorruption',...
        'Number of vertices (%d) does not match expected (%d) !',...
        size(vertices,1),sa/2+2);
end

end
