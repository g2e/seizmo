function takeoffangle=mkrayp2angle(phase,h,p,model,p5,p6);
% mkrayp2ngle........convert ray parameter into take off angle
%
% call: takeoffangle=mkrayp2angle(phase,h,p,model);
%       takeoffangle=mkrayp2angle(phase,h,p,vp,vs,rp);
%
%       phase: string containing seismic phase name, like 'P', 'S', 'ScS', 'PKPdf',...
%              The routione does not test whether your ray parameter makes sense for
%              the phase you desire!
%           h: hypocentral depth [km]
%           p: ray parameter [s/deg]
%       model: model structure as returned by MKREADND
%
%       vp: p wave velocity at focal depth [km/s]
%       vs: s wave velocity at focal depth [km/s]
%       rp: planetary radius [km]
%
%       If H is the depth of a discontinuity, VP and VS may contain two
%       velocity values, the first one for the upper side of the discontinuity
%       and the second one for the lower side.
%       In such cases, the routine will return two angle per input ray parameter!
%
% result: takeoffangle: take off angle at source [deg]
%                       Since P is only the horizontal slowness, not the
%                       slowness vector, it is not possible to distinguish
%                       between upgoing and downgoing rays by the value of
%                       P alone. The take of angle is therefore not
%                       uniquely defined by the ray parameter, but has two
%                       solutions for the two cases.
%                       TAKEOFFANGLE will thus consist of twice as many
%                       elements as P, the first half of the list is for the
%                       upgoing rays, and the second half of the list is
%                       for the downgoing ones.
%                       However, if H==0, an upgoing ray does not make any
%                       sense and is therefore removed from the list. (also
%                       because of compatibility issues with previous
%                       version)
%
%                       Complex take off angles will be set to NaN.
%
% using VP and VS as input parameter, you don't need to call time wasting 
% MKINTERPMODEL here, making everythig quicker.
% 
% This routine is derived fro MKANGLE2RAYP
%
% Martin Knapmeyer, 24.02.2004, 27.06.2006

% 27.06.2006: vectorization: lists of ray parameters are allowed, and two
%             velocity values at discontinuities are allowed.

%%% a useful constant
radian=pi/180;

%%% a magic number: if the imaginary part of an angle is smaller than this,
%%% it is considered to be a numerical artefact and removed.
imageps=1e-4;


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
   %%% velocities are given as input, copy them into the right variables!
   vp=model;
   vs=p5;
   rp=p6;
end; % if nargin

%%% which velocity is relevant?
switch starttype
   case {'s'}
       v=vs;
   case {'p'}
       v=vp;
   otherwise
      error(['MKRAYP2ANGLE: phase names cannot start with ' starttype]);
end; % switch


%%% the computation MK27062006
if length(p)>1
    p=p(:);
end; % if length(p)>1
switch length(v)
    case {1}
         %%% velocity is defined uniquely at depth
         downgoing=asin(p*v/(radian*(rp-h)))/radian;
         upgoing=180-downgoing;
    case {2}
         %%% source is at a velocity jump, upgoing rays see another
         %%% velocity than downgoing rays.
         upgoing=180-asin(p*v(1)/(radian*(rp-h)))/radian;
         downgoing=asin(p*v(2)/(radian*(rp-h)))/radian;
         
    otherwise
         %%% more than 2 samples at the same depth are not allowed!
         error(['MKRAYP2ANGLE: overdetermined velocity at z=' num2str(h)...
                ' in ' focus.name...
                ' (' int2str(length(v)) 'samples).']);
end; % switch length(v)
takeoffangle=[upgoing; downgoing];


%%% find toatally reflected rays for which angle==90
%%% the angle computed above probably has a small imaginary component,
%%% which is to be removed! MK26062006
indies=find((real(takeoffangle)==90)&(abs(imag(takeoffangle))<imageps));
takeoffangle(indies)=real(takeoffangle(indies));


%%% remove the really complex angles by setting them to NaN MK27062006
indies=find(abs(imag(takeoffangle))>imageps);
takeoffangle(indies)=takeoffangle(indies)*0+NaN;

%%% remove upgoing rays from solution for surface foci MK27062006
if h==0
   takeoffangle=takeoffangle((length(p)+1):end);
end; % if h==0

%%% everything that is returned is real, after weh have removed the
%%% imaginary parts. This line just changes the data type.
takeoffangle=real(takeoffangle);


return;