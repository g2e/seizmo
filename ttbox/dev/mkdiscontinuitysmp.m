function [newdepths,criticalz]=mkdiscontinuitysmp(model,zepsilon)
% mkdiscontinuitysmp.......construct new depth samples based on discontinuities
%
% call: [newdepths,criticalz]=mkdiscontinuitysmp(model);
%       [newdepths,criticalz]=mkdiscontinuitysmp(model,zepsilon);
%
%       model: MODEL structure as returned by MKREADND or MKCLR2MODEL
%
%    zepsilon: distance from important points in model to add new points
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
% Garrett Euler,    18.02.2008

%%% a magical number (in km)
%%% vertical epsilon value for placing new samples around important parts
%%% of the velocity model
if(nargin==1); zepsilon=0.001; end

%%%%%% first step: identify discontinuities in OLDMODEL
%%%%%%             As discontinuity, we define any case in which the same depth sample
%%%%%%             appears twice. This might be the case at velocity jumps or at places
%%%%%%             where v(z) has a kink and the first derivative jumps.
%%%%%%             This is not done by MKREADND, so we have to do it here.
%%%%%%             The identification is simple: whenever the difference between
%%%%%%             adjacent depth samples is zero, we have a discontinuity.
%%%%%%             We assume that all depth samples are ordered by depth
%%%%%%             properly!
idiscon=(~diff(model.z)); % upper part of discon
criticalz=model.z([idiscon; false]);

%%%%%% second step: add new samples some epsilon above and below discontinuities 
%%%%%%              of velocity (here we define only where to add samples,
%%%%%%              an interpolation of the model to ALL new samples is done later!)
newdepths=[criticalz-zepsilon; criticalz+zepsilon];

end
