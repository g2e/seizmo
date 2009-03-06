function [newdepths,criticalz]=mkdiscontinuitysmp(model);
% mkdiscontinuitysmp.......construct new depth samples based on discontinuities
%
% call: [newdepths,criticalz]=mkdiscontinuitysmp(model);
%
%       model: MODEL structure as returned by MKREADND or MKCLR2MODEL
%
% result: newdepths: depths of the depth samples to be added to the model [km]
%                    For each discontinuity, more than one (namely two) new
%                    depth samples are generated.
%
%         criticalz: list of the critical depths [km]
%                    This list contains all the critical depths: the actual
%                    depth samples of the discontinuities themselves.
%                    Critical rays can must derived from these.
%
% The flat earth transform distorts the velocity gradient at bot sides of
% discontinuities. In order to predicht the critical ray parameters that
% mark the beginning and and of the total reflection zone, it is necessary
% to repair this. A sufficient method is to insert one sample slightly
% above and another one slightly below the discontinuity.
%
% The new depth samples generated here are intended to be inserted by a
% call of MKINTERPMODEL.
%
% Martin Knapmeyer, 23.06.2003


%%%%%% zeroth step: some initializations
criticalz=[]; % to collect critical depths detected in spherical co'
newdepths=[]; % to collect all extra samples


%%% a magical number
%%% vertical epsilon value for placing new samples around important parts
%%% of the velocity model
zepsilon=0.001; % [km]

%%%%%% first step: identify discontinuities in OLDMODEL
%%%%%%             As discontinuity, we define any case in which the same depth sample
%%%%%%             appears twice. This might be the case at velocity jumps or at places
%%%%%%             where v(z) has a kink and the first derivative jumps.
%%%%%%             This is not done by MKREADND, so we have to do it here.
%%%%%%             The identification is simple: whenever the difference between
%%%%%%             adjacent depth samples is zero, we have a discontinuity.
%%%%%%             we assume that all depth samples are ordered by depth
%%%%%%             properly!
z=model.z;
deltaz=diff(z);
disconindies=find(deltaz==0); % array index to discontinuities: points to the upper sample
disconcnt=length(disconindies); % number of discontinuities
discondepths=z(disconindies); % depths of discontinuities
criticalz=[criticalz; z([disconindies; disconindies+1])]; % to collect all critical depths


%%%%%% second step: add new samples some epsilon above and below first order
%%%%%%              discontinuities of velocity (here we define only where to add samples,
%%%%%%              an interpolation of the model to ALL new samples is done later!)
%%% distance between new samples and the discontinuity
%%% construct extended depth samples list
disconextension=[discondepths-zepsilon; discondepths+zepsilon];
%disp(['MKDISCONTINUITYSMP: ' int2str(length(disconindies))...
%      ' critical depths from discontinuities identified.' ]);
  

%%%%%% third step: combine all new samples identified in spherical co'
criticalz=criticalz;
newdepths=[disconextension];