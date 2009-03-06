function [newdepths,criticalz]=mkderivsmp(oldmodel,zepsilon,thresh)
%mkderivsmp    analyze first derivative of v(z) to find triplications
%
% call: [newdepths,criticalz]=mkderivsmp(oldmodel);
%       [newdepths,criticalz]=mkderivsmp(oldmodel,zepsilon);
%       [newdepths,criticalz]=mkderivsmp(oldmodel,zepsilon,thresh);
%
%       model: MODEL structure as returned by MKREADND or MKCLR2MODEL
%
% result: newdepths: depths of the depth samples to be added to the model.
%                    For each discontinuity, more than one (namely two) new
%                    depth samples are added. CURRENTLY NOT IMPLEMENTED.
%
%
%         criticalz: list of the critical depths [km]
%                    This list contains all the critical depths: the actual
%                    depth samples of the discontinuities themselves.
%                    Critical rays must be derived from these.
%
% Strong changes in the velocity gradient produce triplications in the
% travel time. The beginning (=highest ray parameter) of these
% triplications can be derived from the velocity model, whereas their
% endings (=smallest ray parameter) depend on seismic phase and focal
% depth.
%
% The analysis has to be done in flat earth coordinates.
%
% The new depths (possibly) generated here are intended for insertion into 
% the velocity sampling by a call of MKINTERPMODEL.
%
% Martin Knapmeyer, 23.06.2003
% Garrett Euler, 18.02.2008

%%% magical numbers
if(nargin<2); zepsilon=0.001; end % (km) new depth values distance from critical depths
if(nargin<3); thresh=1.04; end % default minimum ratio of gradients expected to trigger triplication

warning off MatLab:divideByZero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               FLAT EARTH CORDINATES SECTION                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The analysis should be done in spherical coordinates and another part
%%% in flat earth coordinates. using the COORDSYS flag, the flat earth
%%% transformation can be switched off for testing.  This option was 
%%% introduced for testing purposes, but remained in the code.
coordsys='flat';


%%%%%% FET: transformation into flat earth (can be deactivated) since TTBOX
%%%%%%      does its ray tracing in flat earth, the properties of the flat
%%%%%%      earth velocity distribution have to be analyzed!
switch coordsys
    case {'flat'}
        [z,v]=mksfer2flat(oldmodel.rp,oldmodel.z,[oldmodel.vp oldmodel.vs]);
        mkflat2sfer(oldmodel.rp,z)
    case {'spherical'}
        v=[oldmodel.vp oldmodel.vs];
        z=oldmodel.z;
    otherwise
        error('flat or spherical only')
end;


%%%%%% derivatives
dz=diff(z);
dvpdz=diff(v(:,1))./dz;
dvsdz=diff(v(:,2))./dz;


%%%%%% third step:  Identifies strong velocity gradient changes and
%%%%%%              discontinuities.
vpstronggrad=((dvpdz(1:end-1)./dvpdz(2:end))>thresh | (dvpdz(1:end-1)./dvpdz(2:end))<1/thresh);
vsstronggrad=((dvsdz(1:end-1)./dvsdz(2:end))>thresh | (dvsdz(1:end-1)./dvsdz(2:end))<1/thresh);


%%%%%% final step:  Spherical depths of strong gradient changes
criticalz=oldmodel.z([false; vpstronggrad | vsstronggrad; false]);
newdepths=[criticalz-zepsilon; criticalz+zepsilon];


warning on MatLab:divideByZero
end
