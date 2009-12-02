function [a,b,c]=sph_tri_auto(maxdeg)
%SPH_TRI_AUTO    Automatically create spherical triangle polyhedra
%
%    Usage:    [a,b,c]=sph_tri_auto(maxdeg)
%
%    Description: [A,B,C]=SPH_TRI_AUTO(MAXDEG) automatically creates a
%     spherical triangle polyhedra where the vertices are all approximately
%     less than MAXDEG degrees apart.  Distortion due to projection of the
%     triangle vertices onto the unit sphere may lead to some vertices
%     exceeding MAXDEG.
%
%    Notes:
%     - Starts with an icosahedron via SPH_TRI_INIT
%     - The nearest vertices on an icosahedron are 63.435 deg
%
%    Examples:
%     Make a polyhedra for Earth requiring ~50km spacing between vertices
%     (warning this makes a >400000 faceted polygon!):
%      [a,b,c]=sph_tri_auto(50*180/(6371*pi));
%      a=a*6371; b=b*6371; c=c*6371;
%
%     Polyhedra with ~4deg vertex spacing:
%      [a,b,c]=sph_tri_auto(4);
%
%    See also: SPH_TRI_INIT, SPH_TRI_SPLIT, SPH_TRI_VERTICES,
%              SPH_TRI_LATLON, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov. 12, 2009 - initial version
%        Nov. 17, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 17, 2009 at 07:45 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

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
[a,b,c]=sph_tri_init;
[a,b,c]=sph_tri_split(a,b,c,numel(s),s);

end
