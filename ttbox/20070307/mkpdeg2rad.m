function prad=mkpdeg2rad(pdeg);
% mkpdeg2rad...........transform ray parameter from s/deg to s/rad
%
% call: prad=mkpdeg2rad(pdeg);
%
%       pdeg: ray parameter in sec/degree, as usually given in
%             traveltime tables
%
% result: prad: ray parameter in sec/radian, as needed for ray computations
%
% the transformation is done by mutiplying pdeg with 180/pi.
% This function is written simply to  make clear that the transformation is
% important and to show how it works.
%
% Martin Knapmeyer, 22.04.2002


prad=180*pdeg/pi;