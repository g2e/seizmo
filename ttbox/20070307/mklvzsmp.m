function [lvzextra,criticalz]=mklvzsmp(oldmodel);
% mklvzsmp..........analyze model for low velocity zones and their effects
%
% call: [lvzextra,criticalz]=mklvzsmp(oldmodel);
%
%
%
% result: lvzextra: collection of new depth samples derived from LVZ
%                   analysis.
%                   This has the form of a model structure, but contains
%                   only the vp, vs, rho, qp, qs, z values of newly
%                   constructed depth samples.
%                   To expand a model structure coherently, is is necessary
%                   to expand all five model paramaters!
%         criticalz: list of the critical depths [km]
%                    This list contains all the critical depths: the actual
%                    depth values of the LVZ bottoms.
%                    Critical rays must be derived from these.
%
%
% Low velocioty zones cause shadow zones and possibly other effects,
% depending on focal depth and details of the LVZ structure. This routine
% identifies LVZ in the model and returns critical depths and possible
% additional depth samples to improve the model sampling.
%
% Martin Knapmeyer, 20.07.2006, 25.09.2006, 27.09.2006


%%% init result
lvzextra=mkemptymodel;
criticalz=[];


%%% get new samples and critical depths for Vp
[newvp,newvpsec,newzp,critzp]=mklvzanalysis(oldmodel.vp,oldmodel.vs,oldmodel.z,oldmodel.rp);


%%% get new samples and critical depths for Vs
[newvs,newvssec,newzs,critzs]=mklvzanalysis(oldmodel.vs,oldmodel.vp,oldmodel.z,oldmodel.rp);


%%% join critical depths
criticalz=[critzp; critzs];


%%% interpolate model at NEWZP and NEWZS depths to obtain rho and Q there
lvzextra=mkinterpmodel(oldmodel,unique([newzp; newzs]),'simple');



%%% overwrite interpolated velocities with those obtained from LVZ analysis
%%% keep only z, rho and Q values.
%%% first copy values from analysis of Vp profile
newzplen=length(newzp);
for newzpcnt=1:newzplen
    indy=find(lvzextra.z==newzp(newzpcnt));
    lvzextra.vp(indy)=newvp(newzpcnt);
    lvzextra.vs(indy)=newvpsec(newzpcnt);
end; % for newzpcnt
%%% the copy values from analysis of Vs profile
newzslen=length(newzs);
for newzscnt=1:newzslen
    indy=find(lvzextra.z==newzs(newzscnt));
    lvzextra.vp(indy)=newvssec(newzscnt);
    lvzextra.vs(indy)=newvs(newzscnt);
end; % for newzpcnt
%%% everything copied.






