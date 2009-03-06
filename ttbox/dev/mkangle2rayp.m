function p=mkangle2rayp(phase,h,angle,model,p5,p6);
% mkangle2rayp........convert take off angle into ray parameter
%
% call: p=mkangle2rayp(phase,h,angle,model);
%       p=mkangle2rayp(phase,h,angle,vp,vs,rp);
%
%       phase: string containing seismic phase name, like 'P', 'S', 'ScS', 'PKPdf',...
%              The routione does not test whether your ray partameter makes sense for
%              the phase tyou desire!
%           h: hypocentral depth [km]
%       angle: take off angle at source [deg]
%       model: model structure as returned by MKREADND
%
%       vp: p wave velocity at focal depth [km/s]
%       vs: s wave velocity at focal depth [km/s]
%       rp: planetary radius [km]
%
%       VP and VS have to scalar values, the discontinuity-decision depending
%       on the angle has to be made already!
%
% result: p: ray parameter [s/deg]
%
% using VP and VS as input parameter, you don't need to call time wasting MKINTERPMODEL
% here, making everythig quicker.
%
% Martin Knapmeyer, 25.09.2003, 14.11.2006

% 14.11.2006: vectorized velocity choice, depending on takeoff angle: if at
%             a discontinuity, the velocity for upgoing and downgoing rays
%             must be different!

%%% a useful constant
radian=pi/180;

%%% what kind of wave do we emit?
starttype=lower(phase(1));

%%% determine velocities from model structure, if necessary
if nargin==4
   %%% determine velocity at source depth (needed for angle-ray parm conversion)
   %%% if h is a discontinuity depth, there will be two velocioty values.
   %%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
   focus=mkinterpmodel(model,h,'simple');
   indy=find(focus.z==h);
   vp=focus.vp(indy);
   vs=focus.vs(indy);
   rp=model.rp;
else
   %%% velocities are given by the user, just redistribute input over
   %%% the respective variables
   vp=model;
   vs=p5;
   rp=p6;
end; % if nargin


%%% which velocity is relevant? This depends only on first ray leg.
switch starttype
   case {'s'}
       v=vs;
   case {'p'}
       v=vp;
   otherwise
      error(['MKANGLE2RAYP: phase names cannot start with ' starttype]);
end; % switch


%%% at a discontinuity, there will be two velocity values.
%%% use above-value for angle>=90 and below-value for angle<90
%%% it is assumed that the first value is the above-value.
if length(v)==2
  %%% we're at a discontinuity: distinguish between upgoing and
  %%% downgoing rays
  whichangle=(angle<=90)+1;
else
  %%% not at a discontinuity: always use first velocity
  whichangle=zeros(size(angle))+1;
end; % if length(vp)==2


%%% construct v vector according to whichangle list, if necessary
if length(v)==2
    vlist=whichangle;
    vlist=v(whichangle);
else
    vlist=v;
end; % if length(v)==2


%%% the computation
p=radian*sin(angle*radian)*(rp-h)./vlist;


return;