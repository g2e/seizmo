function pdeg=mkprad2deg(prad);
% mkprad2deg...........transform ray parameter from s/rad to s/deg
%
% call: pdeg=mkprad2deg(prad);
%
%       prad: ray parameter in sec/radian, as needed for ray computations
%
% result: pdeg: ray parameter in sec/degree, as usually given in
%               traveltime tables
%
% the transformation is done by mutiplying pdeg with 180/pi.
% This function is written simply to  make clear that the transformation is
% important and to show how it works.
%
% Martin Knapmeyer, 22.04.2002


pdeg=pi*prad/180;