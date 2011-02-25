function [v,tri]=sph_tri_init(polygon)
%SPH_TRI_INIT    Returns vertices for regular spherical triangle polyhedra
%
%    Usage:    [vertices,tri]=sph_tri_init()
%              [vertices,tri]=sph_tri_init(polygon)
%
%    Description:
%     [VERTICES,TRI]=SPH_TRI_INIT() returns the vertices of an icosahedron
%     in VERTICES and the vertex indices for each face in TRI.  VERTICES
%     is a 12x3 array with the columns cooresponding to x y z coordinates
%     and each row belongs to a separate vertex of the icosahedron.  TRI
%     is a 20x3 array of indices going from 1 to 12 indicating the triple
%     of vertices that make each of the 20 faces.  Coordinates are given
%     such that all vertices lie on the unit circle.
%
%     [VERTICES,TRI]=SPH_TRI_INIT(POLYGON) allows specifying an alternative
%     polygon to return.  POLYGON must be 'tetrahedron', 'octahedron' or
%     'icosahedron'.  The default is 'icosahedron'.
%
%    Notes:
%
%    Examples:
%     % Returns the vertices associated with the indicated polygons:
%     sph_tri_init('tetrahedron')
%     sph_tri_init('octahedron')
%     sph_tri_init('icosahedron')
%
%     % Plot an icosahedron:
%     [v,tri]=sph_tri_init;
%     tri=tri(:,[1:3 1]); % close triangles
%     plot3(v(tri)',v(12+tri)',v(24+tri)');
%     axis vis3d
%     axis off
%     rotate3d on
%
%    See also: SPH_TRI_SPLIT, SPH_TRI_AUTO, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov.  6, 2009 - initial version
%        Nov. 16, 2009 - all triangles proceed clockwise
%        Feb. 11, 2011 - minor doc formatting
%        Feb. 18, 2011 - major revision: use vertex & index output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2011 at 15:05 GMT

% todo:

if(~nargin || isempty(polygon)); polygon='icosahedron'; end
if(~ischar(polygon))
    error('seizmo:sph_tri_init:badPoly',...
        ['POLYGON must be ''tetrahedron'', '...
        '''octahedron'', ''icosahedron''!']);
end
switch lower(polygon)
    case 'tetrahedron'
        p=1/sqrt(3);
        v=[p p p; -p -p p; -p p -p; p -p -p];
        tri=[1 2 3; 4 3 2; 3 4 1; 2 1 4];
    case 'octahedron'
        v=[1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1; 0 0 -1];
        tri=[1 5 2; 2 5 3; 3 5 4; 5 1 4; 2 6 1; 3 6 2; 4 6 3; 6 4 1];
    case 'icosahedron'
        p=2*cos(pi/5);
        v=[0 p 1;  0 -p 1; 0  p -1;  0 -p -1; ...
           1 0 p; -1  0 p; 1  0 -p; -1  0 -p; ...
           p 1 0; -p  1 0; p -1  0; -p -1  0];
        nv=sum(abs(v).^2,2).^(0.5);
        v=v./nv(:,[1 1 1]);
        tri=[2 11  4;  5 11 2;  9 11 5;  7 11 9; 11  7 4; 4 12 2; ...
             6  2 12;  2  6 5;  1  5 6;  5  1 9;  3  9 1; 9  3 7; ...
             8  7  3;  7  8 4; 12  4 8; 12 10 6;  6 10 1; 1 10 3; ...
             3 10  8; 10 12 8];
    otherwise
        error('seizmo:sph_tri_init:badPoly',...
            ['POLYGON must be ''tetrahedron'', '...
            '''octahedron'', ''icosahedron''!']);
end

end
