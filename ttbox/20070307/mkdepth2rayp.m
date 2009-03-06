function [prayp,srayp,newz]=mkdepth2rayp(model,z);
% mkdepth2rayp........compute rayp parameters of rays turning at given depth
%
% call: [prayp,srayp,newz]=mkdepth2rayp(model,z);
%
%       model: MODEL structure containing velocity model
%              Thsi is a structure as returned by MKREADND.
%           z: list of depths for which ray parameters are to be computed
%
% result: prayp: ray parameters of P wave rays that have vertex at depths Z
%         srayp: as PRAYP, but for S waves
%         newz:  depth corrsposnding to the ray parameter lists [km]
%                The lengths of PRAYP and SRAYP may be different of that of
%                Z, therefore we return a fitting list here.
%
% A ray with ray parameter as given in PRAYP and SRAYP will travel horizontally
% at depth given in Z. prayp(i) and srayp(i) correspond to z(i).
%
% Martin Knapmeyer, 06.11.2003


%%% a useful constant
radian=pi/180;

%%% interpolate model to depths given in Z
newmodel=mkinterpmodel(model,z,'simple');


%%% the computation
%%% p=radian*sin(angle*radian)*(rp-h)/v;
%%% where angle=90 and v is taken from NEWMODEL.
%%% note that v might be zero, so we have to take care to avoid divide by zero
zeroelements=find(newmodel.vp==0);
nonzeros=find(newmodel.vp~=0);
prayp=zeros(size(newmodel.z));
prayp(nonzeros)=radian*(newmodel.rp-newmodel.z(nonzeros))./newmodel.vp(nonzeros);
prayp(zeroelements)=prayp(zeroelements)+NaN;

zeroelements=find(newmodel.vs==0);
nonzeros=find(newmodel.vs~=0);
srayp=zeros(size(newmodel.z));
srayp(nonzeros)=radian*(newmodel.rp-newmodel.z(nonzeros))./newmodel.vs(nonzeros);
srayp(zeroelements)=srayp(zeroelements)+NaN;


%%% construct corresponding depth list
newz=newmodel.z;