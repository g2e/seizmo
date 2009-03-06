function radii=mkdetectdiscon(model);
% mkdetectdiscon.........detect velocity discontinuities from z sample list
%
% call: radii=mkdetectdiscon(model);
%
%       depthlist: velocity model structure as returned by MKREADND
%
% result: radii: radii of all discontinuities
%                empty if no discontinuities exist.
%
% The routine assumes that discontinuities can be identified by the
% existence of double depth samples: if two adjacent samples have the same
% depth, the existence of a discontinuity is concluded.
%
% The result is returned as radius rather than depth since MKRAYDEPTH, for
% which we do this, thinks in terms of radii.
%
% Martin Knapmeyer, 04.12.2006


%%% init result;
radii=[];


%%% find adjacent list elements with identical depth
zdiff=diff(model.z);
zindy=find(zdiff==0);


%%% construct radius list from index list
radii=model.rp-model.z(zindy);