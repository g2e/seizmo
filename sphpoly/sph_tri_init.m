function [varargout]=sph_tri_init(polygon)
%SPH_TRI_INIT    Returns vertices for regular spherical triangle polyhedra
%
%    Usage:    [a,b,c]=sph_tri_init()
%              [a,b,c]=sph_tri_init(polygon)
%              [vertices]=sph_tri_init(...)
%
%    Description:
%     [A,B,C]=SPH_TRI_INIT() returns the vertices of the faces of a
%     spherical icosahedron in A, B, & C.  A, B, & C are 20x3 arrays with
%     the columns cooresponding to x y z coordinates and each row belongs
%     to a separate face of the icosahedron.  So A(1,:), B(1,:), & C(1,:)
%     give the vertices of the first face of the icosahedron.  Coordinates
%     are given such that all vertices lie on the unit circle.
%
%     [A,B,C]=SPH_TRI_INIT(POLYGON) allows specifying an alternative
%     polygon to return vertices for.  POLYGON must be 'tetrahedron',
%     'octahedron' or 'icosahedron'.  The default is 'icosahedron'.
%
%     [VERTICES]=SPH_TRI_INIT(...) returns the unique vertices associated
%     with the polygon rather than replicating them for each face that
%     shares the vertices.  This reduces the output to an Nx3 array where
%     N is NUMBER_OF_FACES/2 + 2.
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
%     [a,b,c]=sph_tri_init;
%     x=[a(:,1) b(:,1) c(:,1) a(:,1)]';
%     y=[a(:,2) b(:,2) c(:,2) a(:,2)]';
%     z=[a(:,3) b(:,3) c(:,3) a(:,3)]';
%     plot3(x,y,z);
%     axis vis3d
%     axis off
%
%    See also: SPH_TRI_SPLIT, SPH_TRI_AUTO, SPH_TRI_VERTICES,
%              SPH_TRI_LATLON, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov.  6, 2009 - initial version
%        Nov. 16, 2009 - all triangles proceed clockwise
%        Feb. 11, 2011 - minor doc formmatting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

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
        idx=[p p p; -p -p p; -p p -p; p -p -p];
        if(nargout<=1); varargout{1}=idx; return; end
        varargout{1}=idx([1 4 3 2],:);
        varargout{2}=idx([2 3 4 1],:);
        varargout{3}=idx([3 2 1 4],:);
    case 'octahedron'
        idx=[1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1; 0 0 -1];
        if(nargout<=1); varargout{1}=idx; return; end
        varargout{1}=idx([1 2 3 5 2 3 4 6],:);
        varargout{2}=idx([5 5 5 1 6 6 6 4],:);
        varargout{3}=idx([2 3 4 4 1 2 3 1],:);
    case 'icosahedron'
        p=2*cos(pi/5);
        idx=[0 p 1;  0 -p 1; 0  p -1;  0 -p -1; ...
             1 0 p; -1  0 p; 1  0 -p; -1  0 -p; ...
             p 1 0; -p  1 0; p -1  0; -p -1  0];
        nidx=sum(abs(idx).^2,2).^(0.5);
        idx=idx./nidx(:,[1 1 1]);
        if(nargout<=1); varargout{1}=idx; return; end
        varargout{1}=...
            idx([ 2  5  9  7 11  4  6 2 1 5 3 9 8 7 12 12  6  1  3 10],:);
        varargout{2}=...
            idx([11 11 11 11  7 12  2 6 5 1 9 3 7 8  4 10 10 10 10 12],:);
        varargout{3}=...
            idx([ 4  2  5  9  4  2 12 5 6 9 1 7 3 4  8  6  1  3  8  8],:);
    otherwise
        error('seizmo:sph_tri_init:badPoly',...
            ['POLYGON must be ''tetrahedron'', '...
            '''octahedron'', ''icosahedron''!']);
end

end
