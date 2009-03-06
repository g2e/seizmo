function p=mkpoisson(vp,vs);
% mkpoisson.........compute poisson ratio from Vp and Vs
%
% call: p=mkpoisson(vp,vs);
%
%         vp: P wave velocity
%         vs: S wave velocity
%
%         vp and vs have to be scalars
%
% result: p: poisson ratio
%
% p is 0 for vp/vs=sqrt(2)
% p is 0.25 for vp/vs=sqrt(3)
% p is 0 for vs=0 (vp/vs=inf)
%
% Martin Knapmeyer 12.04.2002, 09.12.2002, 15.09.2003

discr=vp.*vp-vs.*vs;
if discr~=0
   p=(0.5.*vp.*vp-vs.*vs)./discr;
else
   p=NaN; % this happens only if vp==vs - hence never with real velocities!
end; % if discr

% if vs==0
%    p=0.5; % vs==0 in fluids, p=0.5 then
% else
%    vpvs2=(vp*vp)/(vs*vs);
%    discr=vpvs2-1;
%    if discr==0
%       p=NaN; % this happens only if vp==vs - hence never with real velocities!   
%    else
%       p=(vpvs2-2)/(2*(vpvs2-1));
%    end; % if discr==0 else
% end; % if vs==0 else