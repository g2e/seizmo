function [lat,lon]=sph_tri_latlon(a,b,c)
%SPH_TRI_LATLON    Convert spherical triangle xyz arrays to lat/lon
%
%    Usage:    [lat,lon]=sph_tri_latlon(a,b,c)
%
%    Description:
%     [LAT,LON]=SPH_TRI_LATLON(A,B,C) converts a spherical triangle
%     polyhedra from xyz form to latitude and longitude.  A, B, & C are the
%     xyz matrices of size Nx3 where N is the number of triangles of the
%     polyhedra.  LAT & LON are 3xN where the first row corresponds to the
%     vertices in A, the 2nd row corresponds to B and the 3rd row gives the
%    lat/lon for C.
%
%    Notes:
%
%    Examples:
%     % Get the areas of the triangles for a polyhedra:
%     [a,b,c]=sph_tri_init;
%     [lat,lon]=sph_tri_latlon(a,b,c);
%     sph_poly_area(lat,lon)'
%
%    See also: SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_TRI_AUTO, SPH_TRI_VERTICES,
%              SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov. 16, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix, minor doc formmatting
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
    error('seizmo:sph_tri_split:badInput',...
        'A, B, & C must be equal-sized Mx3 arrays!');
elseif(any(abs(vecnorm(cat(1,a,b,c),2)-1)>sqrt(3*eps*eps)))
    error('seizmo:sph_tri_split:badInput',...
        'A, B, & C vertices must be on the unit sphere!');
end

% convert to xyz
x=cat(1,a(:,1)',b(:,1)',c(:,1)');
y=cat(1,a(:,2)',b(:,2)',c(:,2)');
z=cat(1,a(:,3)',b(:,3)',c(:,3)');

% convert to lat/lon
[lat,lon]=xyz2geocentric(x,y,z);

end
