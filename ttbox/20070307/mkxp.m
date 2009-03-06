function x=mkxp(p,v,z,zmin,zmax,novertex);
% mkxp.........solve integral for x(p) IN FLAT EARTH
%
% call: x=mkxp(p,v,z);
%       x=mkxp(p,v,z,zmin,zmax);
%       x=mkxp(p,v,z,zmin,zmax,novertex);
%
%         p: ray parameter [sec/km]
%         v: wave velocities at upper and lower boundary of layer [km/s]
%            V(1): upper boundary, V(2): lower boundary
%         z: depths of upper and lower boundary of layer [km]
%            Z(1): upper boundary, Z(2): lower boundary
%         zmin: minimum depth within the layer [km]
%               by using zmin, the case that a source is in the layer is handeled
%               If not given, ZMIN is set to Z(1) internally.
%               If ZMIN<Z(1), it is set to Z(1) internally
%         zmax: maximum depth reached by ray within the layer [km]
%               If not given, it is assumed that the ray will penetrate
%               a deeper layer, thus ZMAX is set to Z(2) internally.
%               If ZMAX>Z(2), ZMAX is set to Z(2) internally, what makes handling of
%               rays transmitted to the next layer easier. (hopefully)
%               If Z(1)<ZMAX<Z(2) it is assumed that ZMAX is the depth of the ray
%               vertex (not tested to save computation time). Then another formula
%               is used that allows precise calculation at the vertex.
%         novertex: if present, no vertex is assumed independent of the value
%               of ZMAX. This is needed for focal depth corrections for
%               non-surface foci. The value of MODE is ignored.
%
%
% result: x: horizontal distance that the ray travels within layer [km]
%            - NaN will be returned if the computed x has an imaginary part - this is the
%              case, when ZMAX is deeper than the vertex depth for P
%            - x=0 if z(1)==z(2) - no computation is done then.
%            - inf will be returned if velocity is zero within the layer.
%
% IMPORTANT: this function computes the horizontal distance betwee ray at depth ZMIN and
% ray at depth ZMAX. If the ray has a vertex within the layer, ZMAX has to be the vertex
% depth and the resulting X has to be doubled to get the total distance that the ray
% travels within the layer!
%            
%
% rays may travel fdrom boundary to boundary or from a source within the layer to a boundary
%
% THIS ROUTINE IS NOT VECTORIZED. ONLY ONE LAYER IS HANDELED PER CALL.
%
% structure of this routine is created analogous to tha of MKTP.
%
% Martin Knapmeyer, 22.04.2002, 24.04.2002, 30.04.2002, 12.06.2002


%disp(['MKXP: starting. p=' num2str(p)]);

%%% check if z(1)==z(2)
%%% this would be the case at discontinuities
%%% if so, return t=0
if z(1)==z(2)
   %disp('MKXP: discontinuity!');
   x=0;
   return;
end; % if z(1)
%%% otherwise continue with calculating


%%% check if velocity is constant zero within layer
%%% if so, return x=inf
if (v(1)==0)&&(v(2)==0)
   %disp('MKXP: zero velocity!');
   x=inf;
   return;
end; % if v(1)
%%% otherwise continue with calculating


%%% check input
nin=nargin;
if nin==3
   zmin=z(1);
   zmax=z(2);
else
   zmax=min(z(2),zmax);
   zmin=max(z(1),zmin);

   if nin==5
      novertex=0; % 0: allow execution of vertex code
   else
      novertex=novertex; % use NOVERTEX as given byr user (OK, this is trivial, but...)
   end; % if nin==5
end; % if nin==3


%%% the computation itself
%%% the integral will be evaluated at points ZMIN and ZMAX
%%% for this evaluations, the parameters g and v0 of the linear velocity law are needed
[g,v0]=mkgv0(v,z);
if g==0
   %%% for media with constant v, a more simple solution is used
   x=-p*v0*(zmin-zmax)/sqrt(1-p*p*v0*v0);
else
   %%% media with velocity gradient
   if (z(1)<zmax)&&(zmax<z(2))&&(novertex==0) % was:(~(novertex==1))
      %%% DO NOT CHANGE THE CONDITIONS ABOVE!!!
      %%% this is the vertex code suppressed by NOVERTEX
      %%% vertex at depth zmax, which is between z(1) and z(2)
      %%% recompute g and v0 for z=z-z(1)
      %%% this is needed since the computation of t below assumes that ray begins at z=0 with v(z)=v0
      %%% g remains the same and v0=v(1)
      v0=v(1);
      %%% helper
      %pv0=p*v0;
      %%% compute distance
      zs=zmax-z(1);
      v0gzs=v0+g*zs;
      x=p*sqrt(g*zs*(v0+v0gzs))*v0gzs/g;
   else
      %%% no vertex
      
      %%% some helpers needed at both ZMIN and ZMAX
      p2g2=p*p*g*g;
      v0rcpg=v0/g;
      %%% evaluation at ZMIN, p2g will be applied later for numerical reasons
      zplusv0g=zmin+v0rcpg;
      xmin=sqrt(-p2g2*zplusv0g*zplusv0g+1);
      %%% evaluation at ZMAX
      zplusv0g=zmax+v0rcpg;
      xmax=sqrt(-p2g2*zplusv0g*zplusv0g+1);
      %%% result of integration
      x=-(xmax-xmin)/(p*g);
          
      if isnan(x)
         error('MKXP: distance computation resulting in NaN!');
      end; % if isnan
      %%% if the vertex is at the lower boundary of the layer, the x computed here will be
      %%% complex. Then we make a new computation with the vertex-code
      %%% asking for imag(x) is not quite correct, since x may be complex e.g. for a PKP that
      %%% misses the core. It would be better to compute the real vertex depth here and compare it
      %%% with z(2).
      %%% This requires to reproduce the vertex-code here. In a recursive call of MKXP,
      %%% the zmax<z(2) condition would avert execution of the vertex code. MK12.06.2002
      if imag(x)~=0
         %%% this is the vertex code from above in compact form
         %disp(['MKXP: recomputing the vertex! p=' num2str(p)]);
         v0=v(1);
         %pv0=p*v0;
         zs=zmax-z(1);
         v0gzs=v0+g*zs;
         x=p*sqrt(g*zs*(v0+v0gzs))*v0gzs/g;         
      end; % if imag(x)
   end; % if (z(1)<zmax)&(zmax<z(2))&(~(novertex==1))
end; % if g==0



%%% handle complex/imaginary results
%%% x will be complex or imaginary if zmax is larger (deeper) than the vertex depth
%%% or if the no-vertex code is evaluated at the vertex.
if imag(x)~=0
   x=NaN;
   %disp(['MKXP: complex distance, returning ' num2str(x) ' instead!']);
end; % if imag(x)

%disp(['MKXP: result: ' num2str(x)]);

%%% return result
% x=x;  % result is already in the output variable
