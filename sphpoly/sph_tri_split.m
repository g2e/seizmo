function [a,b,c]=sph_tri_split(a,b,c,n,s)
%SPH_TRI_SPLIT    Splits up faces of spherical triangle polyhedra
%
%    Usage:    [a,b,c]=sph_tri_split(a,b,c,n)
%              [a,b,c]=sph_tri_split(a,b,c,n,s)
%
%    Description:
%     [A,B,C]=SPH_TRI_SPLIT(A,B,C,N) splits the spherical triangle
%     polyhedra faces defined by A, B, & C (vertex arrays of xyz values -
%     see SPH_TRI_INIT for more details) into 4 approximately equal-area,
%     equalaterial sub-triangles whos vertices are all on a unit sphere.
%     N gives the number of spliting iterations and must be an integer >=0.
%     An N of 2 would split each face of the original polyhedra into 4 and
%     then each of those faces would be further split into 4 sub-triangles.
%
%     [A,B,C]=SPH_TRI_SPLIT(A,B,C,N,S) specifies the splitting order(s) S.
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
%     % Create a regular icosahedron and split each face into 9 triangle on
%     % the unit sphere:
%     [a,b,c]=sph_tri_init;
%     [a,b,c]=sph_tri_split(a,b,c,1,3);
%
%     % Now plot the 180-faceted polygon:
%     x=[a(:,1) b(:,1) c(:,1) a(:,1)]';
%     y=[a(:,2) b(:,2) c(:,2) a(:,2)]';
%     z=[a(:,3) b(:,3) c(:,3) a(:,3)]';
%     plot3(x,y,z);
%     axis vis3d
%     axis off
%
%    See also: SPH_TRI_INIT, SPH_TRI_AUTO, SPH_TRI_VERTICES,
%              SPH_TRI_LATLON, SPH_POLY_AREA, SPH_POLY_IN

%     Version History:
%        Nov.  7, 2009 - initial version
%        Nov. 12, 2009 - added doc
%        Nov. 16, 2009 - all triangles proceed clockwise
%        Feb. 11, 2011 - mass nargchk fix, minor doc formmatting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:
% - really should use Tegmark's method

% check nargin
error(nargchk(4,5,nargin));

% default splitting order
if(nargin==4 || isempty(s)); s=2; end

% check arguments
sa=size(a); sb=size(b); sc=size(c);
if(~isequal(sa,sb,sc) || numel([sa sb sc])>6 ...
        || sa(2)~=3 || sb(2)~=3 || sc(2)~=3)
    error('seizmo:sph_tri_split:badInput',...
        'A, B, & C must be equal-sized Mx3 arrays!');
elseif(any(abs(vecnorm(cat(1,a,b,c),2)-1)>sqrt(3*eps*eps)))
    error('seizmo:sph_tri_split:badInput',...
        'A, B, & C vertices must be on the unit sphere!');
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
for i=1:n
    % preallocate
    nt=size(a,1); % number of triangles before
    ss=s(i)^2; % number of new triangles per old triangle 
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

end
