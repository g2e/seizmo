function [v,tri]=sph_tri_split(v,tri,n,s)
%SPH_TRI_SPLIT    Splits up faces of spherical triangle polyhedra
%
%    Usage:    [v,tri]=sph_tri_split(v,tri,n)
%              [v,tri]=sph_tri_split(v,tri,n,s)
%
%    Description:
%     [V,TRI]=SPH_TRI_SPLIT(V,TRI,N) splits the spherical triangle
%     polyhedra faces defined by V & TRI (array of vertex xyz values &
%     vectex indices for each face) into 4 approximately equal-area,
%     equalaterial sub-triangles whos vertices are all on a unit sphere.
%     N gives the number of spliting iterations and must be an integer >=0.
%     An N of 2 would split each face of the original polyhedra into 4 and
%     then each of those faces would be further split into 4 sub-triangles.
%
%     [V,TRI]=SPH_TRI_SPLIT(V,TRI,N,S) specifies the splitting order(s) S.
%     S must be a scalar or vector of length N.  Elements of S must be
%     integers >=1 and the default value is 2.  Splitting order is the
%     number of equal-length segments that the edge of each triangle is
%     split into when producing the sub-triangles.  The splitting order is
%     also the square root of the number of sub-triangles produced from a
%     single triangle during that iteration.  So N=2, S=[3 4] for an
%     icosahedron input (20 triangular faces) will produce a spherical
%     triangle polyhedra with 2880 faces (20*9*16).
%
%    Notes:
%
%    Examples:
%     % Create a regular icosahedron and split each face
%     % into 9 triangles on the unit sphere:
%     [v,tri]=sph_tri_init;
%     [v,tri]=sph_tri_split(v,tri,1,3);
%
%     % Now plot the 180-faceted polygon:
%     tri=tri(:,[1:3 1]); % close triangles
%     plot3(v(tri)',v(92+tri)',v(184+tri)');
%     axis vis3d
%     axis off
%
%    See also: SPH_TRI_INIT, SPH_TRI_AUTO, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov.  7, 2009 - initial version
%        Nov. 12, 2009 - added doc
%        Nov. 16, 2009 - all triangles proceed clockwise
%        Feb. 11, 2011 - mass nargchk fix, minor doc formmatting
%        Feb. 18, 2011 - major revision: use vertex & index input/output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 18, 2011 at 15:05 GMT

% todo:
% - another function for Tegmark's method
% - algorithm is poor & uses old input/output
%   - get new vertices
%     - each face separately makes sense
%       - unique issues (need sloppy mode)
%   - build tri (need simple function s2tri)

% check nargin
error(nargchk(3,4,nargin));

% default splitting order
if(nargin==3 || isempty(s)); s=2; end

% check arguments
sv=size(v);
st=size(tri);
if(sv(2)~=3 || numel(sv)>2)
    error('seizmo:sph_tri_split:badInput',...
        'V must be a Nx3 array!');
elseif(any(abs(vecnorm(v,2)-1)>sqrt(3)*eps))
    error('seizmo:sph_tri_split:badInput',...
        'Vertices must be on the unit sphere!');
elseif(st(2)~=3 || numel(st)>2 || st(1)/2+2~=sv(1))
    error('seizmo:sph_tri_split:badInput',...
        'TRI must be a Mx3 array where M=2*(SIZE(V,1)-2) !');
elseif(~isreal(n) || ~isscalar(n) || n~=fix(n) || n<0)
    error('seizmo:sph_tri_split:badInput',...
        'N must be a scalar integer >= 0 !');
elseif(~isreal(s) || (~isscalar(s) && (~isvector(s) || (length(s)~=n))) ...
        || any(s~=fix(s)) || any(s<1))
    error('seizmo:sph_tri_split:badInput',...
        'S must be a scalar or N-element vector of integers >= 1 !');
end

% expand scalar s
if(isscalar(s)); s=s(ones(1,n)); end

% loop over number of steps
[a,b,c]=vtri2abc(v,tri);
for i=1:n
    % preallocate
    nt=size(a,1); % number of triangles before
    ss=s(i)^2;    % number of new triangles per old triangle
    an=nan(nt*ss,3); bn=an; cn=an;
    
    % often used
    s1=s(i)-1;
    s2=s(i)-2;
    sp=s(i)+1;
    
    % for calcing vertices along boundary
    u=1:s1;
    d=fliplr(u);
    k=ones(s1,1);
    
    % for calcing vertices within
    t=1:s2; t=t(:); t=t(:,ones(s2,1));
    l1=tril(true(s2));
    l2=flipud(l1);
    v=t(l2)';
    [y,y]=find(l1); y=y';
    z=s(i)-y;
    w=z-v;
    
    % indices for combining
    ia=1;
    ib=sp;
    ic=(sp*sp+sp)/2;
    ip=2:s(i);
    iq=ic-cumsum(1:s1);
    ir=iq-(1:s1);
    iq=fliplr(iq);
    if(~isempty(y))
        im=s(i)+2*y+(1:sum(1:s2)); % breaks for s=2
    else
        im=[];
    end
    h=nan(ic,3);
    
    % for extracting combined to assigned triangles
    % - cheaper to do with a loop
    x1=nan(1,ss);
    x3=x1;
    midx=false(1,ss);
    x1(1:s(i))=1:s(i);
    x3(1:s(i))=(1:s(i))+sp;
    tmp=cumsum(sp:-1:3);
    c1=s(i); c2=1;
    for ii=s(i):-1:2
        ij=ii-1;
        midx(c1+(1:ij))=true;
        x1(c1+(1:(2*ij)))=tmp(c2)+[1:ij 1:ij];
        x3(c1+(1:(2*ij)))=tmp(c2)+[1:ij 1:ij]+ii*[-ones(1,ij) ones(1,ij)];
        c1=c1+2*ij;
        c2=c2+1;
    end
    x2=x1+1;
    
    % force triangles in the interior to progress clockwise
    x1(midx)=x1(midx)+1;
    x2(midx)=x2(midx)-1;
    
    % loop over triangles
    for j=1:nt
        % vertices on the sides
        p=(d([1 1 1],:)'.*a(j.*k,:)+u([1 1 1],:)'.*b(j.*k,:))./s(i);
        q=(d([1 1 1],:)'.*b(j.*k,:)+u([1 1 1],:)'.*c(j.*k,:))./s(i);
        r=(d([1 1 1],:)'.*c(j.*k,:)+u([1 1 1],:)'.*a(j.*k,:))./s(i);
        
        % project onto unit circle
        pn=vecnorm(p,2);
        qn=vecnorm(q,2);
        rn=vecnorm(r,2);
        p=p./pn(:,[1 1 1]);
        q=q./qn(:,[1 1 1]);
        r=r./rn(:,[1 1 1]);
        
        % combine into a single matrix
        h(ia,:)=a(j,:);
        h(ib,:)=b(j,:);
        h(ic,:)=c(j,:);
        h(ip,:)=p;
        h(iq,:)=q;
        h(ir,:)=r;
        
        % now repeat for interior vertices
        if(~isempty(im))
            m=(w([1 1 1],:)'.*r(z,:)+v([1 1 1],:)'.*q(y,:))./z([1 1 1],:)';
            mn=vecnorm(m,2);
            m=m./mn(:,[1 1 1]);
            h(im,:)=m;
        end
        
        % assign vertices to new triangles
        an(ss*(j-1)+(1:ss),:)=h(x1,:);
        bn(ss*(j-1)+(1:ss),:)=h(x2,:);
        cn(ss*(j-1)+(1:ss),:)=h(x3,:);
    end
    % clear original triangles
    a=an; b=bn; c=cn;
end

% output
[v,tri]=abc2vtri(a,b,c);

end


function [a,b,c]=vtri2abc(v,tri)
a=v(tri(:,1),:);
b=v(tri(:,2),:);
c=v(tri(:,3),:);
end


function [v,tri]=abc2vtri(a,b,c)
v=unique([a;b;c],'rows');
if(size(v,1)~=size(a,1)/2+2)
    warning('seizmo:sph_tri_split:verticesCorruption',...
        'Number of vertices (%d) does not match expected (%d) !',...
        size(v,1),size(a,1)/2+2);
end
[tmp,tri(:,1)]=ismember(a,v,'rows');
[tmp,tri(:,2)]=ismember(b,v,'rows');
[tmp,tri(:,3)]=ismember(c,v,'rows');
end


%function [v]=s2v(v,s)
% from:
%     3
%    1 2
% to (s=2):
%     6
%    4 5
%   1 2 3
% to (s=3):
%     0
%    8 9
%   5 6 7
%  1 2 3 4
%end


%function [tri]=s2tri(s)
%nt=s^2;
%nv=nt/2+2;
% can we break this into simple steps?
%end


