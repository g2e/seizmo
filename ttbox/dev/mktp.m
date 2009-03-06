function t=mktp(p,v,z,zmin,zmax,novertex);
% mktp.........solve integral for traveltime T(p) IN FLAT EARTH
%
% call: t=mktp(p,v,z);
%       t=mktp(p,v,z,zmin,zmax);
%       t=mktp(p,v,z,zmin,zmax,novertex);
%
%         p: ray parameter [sec/km]
%         v: wave velocities at upper and lower boundary of layer [km/s]
%            V(1): upper boundary, V(2): lower boundary
%         z: depths of upper and lower boundary of layer [km]
%            Z(1): upper boundary, Z(2): lower boundary
%         zmin: minimum depth within the layer [km]
%               by using zmin, the case that a source is in the layer is handeled
%               If not given, ZMIN is set to Z(1) internally.
%               If ZMIN<Z(1), it is set to Z(1) internally.
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
% result: t: traveltime of ray between depths ZMIN and ZMAX.
%            set to NaN - something's wrong with your P, ZMIN, or ZMAX then.
%            t=0 if Z(1)==Z(2) - no computation is done then.
%            inf will be returned if velocity is constant zero within layer.
%            inf will also be returned if the pole of the time-integral is hit.
%
% IMPORTANT: this function computes the traveltime between ray at depth ZMIN and
% ray at depth ZMAX. If the ray has a vertex within the layer, ZMAX has to be the vertex
% depth and the resulting T has to be doubled to get the total traveltime that the ray
% spends within the layer!
%
%
% rays may travel from boundary to boundary or from a source within the layer to a boundary
%
% THIS ROUTINE IS NOT VECTORIZED. ONLY ONE LAYER IS HANDELED PER CALL.
%
% Martin Knapmeyer, 24.04.2002, 30.04.2002, 12.06.2002, 19.05.2004


%%% check if z(1)==z(2)
%%% this would be the case at discontinuities
%%% if so, return t=0
if z(1)==z(2)
   t=0;
   return;
end; % if z(1)
%%% otherwise continue with calculating


%%% check if velocity is constant zero within layer
%%% if so, return x=NaN
if (v(1)==0)&(v(2)==0)
   t=inf;
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
   if nin==6 % set flags to suppres or execute vertex code
      novertex=1; % suppress vertex code
   else
      novertex=0; % execute vertex code
   end; % if nin==6
end; % if nin==3




%%% the computation itself
%%% the integral will be evaluated at points zmin and zmax
%%% for this evaluations, the parameters g and v0 of the linear velocity law are needed
[g,v0]=mkgv0(v,z);
if g==0
   %%% for media with constant velocity, the traveltime is simply
   t=(zmax-zmin)/sqrt(v0*(1-p*p*v0)); % not yet tested! mk29.04.2002
else
   %%% media with velocity gradients
   if (z(1)<zmax)&(zmax<z(2))&(~(novertex==1))
      %%% this is the vertex code suppressed by NOVERTEX
      %%% vertex at depth zmax
      %%% recompute g and v0 for z=z-z(1)
      %%% this is needed since the computation of t below assumes that ray begins at z=0 with v(z)=v0
      %%% g remains the same and v0=v(1)
      v0=v(1);
      %%% helper
      pv0=p*v0;
      %%% compute traveltime
      t=log((1+sqrt(1-pv0*pv0))/(pv0))/g;
      %disp(['MKTP: time to vertex ' num2str(t) 's']);
   else
      %%% no vertex
      %%% some helpers needed at both ZMIN and ZMAX
      p2g=p*p*g;
      v0rcpg=v0/g;
      %%% evaluation at ZMIN
      tmin=mkevalt(zmin,g,p2g,v0rcpg);
      %%% evalutation at ZMAX
      tmax=mkevalt(zmax,g,p2g,v0rcpg);
      %%% result of integration
      t=(-1/g)*(tmax-tmin);
      %%% t may be NaN if something's wrong...
      if isnan(t)
         %error('MKTP: time computation resulting in NaN!');
         % NaN is a legal output of this routine, anyybody who calls it must be able to handle NaN!
      end; % if isnan t
      if isinf(real(t))
         t=inf; % just to make sure that it is not inf+iNaN!
      end; % if isinf
      %%% t may be complex here if the vertex is identical with the lower boudnary of
      %%% the layer. Then we make a new computation with the vertex-code. MK 12.06.2002
      if imag(t)~=0
         %%% this is the vertex code from above in compact form
         %disp('MKTP: recomputing the vertex');
         v0=v(1);
         pv0=p*v0;
         t=log((1+sqrt(1-pv0*pv0))/(pv0))/g;
      end; % if imag(t)
   end; % if (z(1)<zmax)&(zmax<z(2) else
end; % if g==0 else


%%% handle complex/imaginary results
%%% t will be complex or imaginary if zmax is larger (deeper) than the vertex depth
if imag(t)~=0
   t=NaN;
end; % if imag(x)

%%% return result
% t=t;  % result is already in the output variable

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                 HELPER FUNCTIONS                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function t=mkevalt(z,g,p2g,v0rcpg);
% mkevalt.......performs the integral evaluation
%
% call: t=mkeavtl(z,g,p2g,v0rcpg);
%
%         z: depth at which the integral is to be evaluated
%         g: velocity gradient as computed by MKGV0
%       p2g: sqr(p)*g, where p is ray parameter and g velocity gradient
%    v0rcpg: v0/g, where v0 is baseline velocity and g velocity gradient (outpout of MKGV0)
%
% result: t: result of traveltime integral at depth z
%
% Traveltime from z1 to z2 is -1/g*(t2-t1)
% the multiplication with 1/g is not applied here for numerical reasons.
%
% Martin Knapmeyer, 26.04.2002

zplusv0g=z+v0rcpg;

denominator=sqrt(-p2g*g*zplusv0g*zplusv0g+1);
if denominator==0 %sqrt(-p2g*g*zplusv0g*zplusv0g+1)==0
   t=0; % this patches a singularity which produced problems in earlier versions. MK19052004
   return;
end;
t=atanh(1/denominator); %sqrt(-p2g*g*zplusv0g*zplusv0g+1));



return;


