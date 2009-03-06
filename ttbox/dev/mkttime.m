function [tt,p]=mkttime(phase,delta,h,model,pscan);
% mkttime.......compute seismic traveltimes
%
% call: [tt,p]=mkttime(phase,delta,h,model);
%       [tt,p]=mkttime(phase,delta,h,model,pscan);
%
%          phase: string containing seismic phase name like 'P', 'S', 'ScS', 'PKPdf', etc.\
%                 Phase names are case sensitive!
%          delta: epicentral distance [deg]
%              h: focal depth [km]
%          model: A structure describing the velocity distribution.
%                 The structure is expected to have the following fields:
%                 model.z: depth [km below surface] (array)
%                 model.vp: Vp [km/s] (array)
%                 model.vs: Vs [km/s] (array)
%                 model.rho: rho [g/ccm] (array)
%                 model.qp: Qp (array)
%                 model.qs: Qs (array)
%                 model.conr: depth of conrad discontinuity
%                 model.moho: depth of moho
%                 model.d410: depth of Mantle-Transition Zone-discontinuity (the "410" on earth)
%                 model.d520: depth of olivine beta-gamma transition (the "520" on earth)
%                 model.d660: depth of lower mantle discontinuioty (the "660" on earth)
%                 model.cmb: depth of core mantle boundary
%                 model.icb: depth of inner core boundary
%                 model.dz: depths of additional discontinuities (array of numbers)
%                 model.dname: names of additional discontinuities (array of strings)
%                 model.rp: planetary radius
%                 model.name: name of model (string)
%                 such a structure can be obtained via MKREADND.
%
%                 if MODEL does not contain the .criticalrays field
%                 created by MKIMPROVEMODEL, MKIMPROVEMODEL is called to
%                 create one. This field is necessary for MKSHOOTRAY.
%
%          pscan: optional argument: a PSCAN structure as returned by
%                 MKPSAMPLER. This may be used to speed up the evaluation.
%
%
% result: tt: traveltime of phase PHASE from a source at distance DELTA and depth H.
%             Due to triplications and the like, this may be a vector!
%             NaN if no time could be computed.
%          p: ray parameters that corrspond to the times in TT 
%             (may be more than one due to triplications at discontinuities!)
%             p(1) yields TT(1) and so on.
%             NaN if MKFINDP could not find a ray parameter
%
% Martin Knapmeyer, 08.05.2002, 05.07.2002, 29.09.2003, 26.08.2004,
%                   19.12.2006

%%% 19.12.2006  use of MKSHOOTRAY implemented


%%% does the IMODEL structure already contain the critical ray paramaters
%%% list? If not, create one! MK1912006
if ~isfield(model,'criticalrays')
    %%% imodel does not contain a cirtical rays list, create one!
    model=mkimprovemodel(model);
end; % if ~isfield(model,'criticalrays')


%%% find ray parameter at which phase arrives at distance DELTA.
%%% p may contain more than one value!
%%% calls mkshootray now MK19122006
if nargin==6
    [p,angles,dists]=mkshootray(phase,delta,h,model,pscan);
else
    [p,angles,dists]=mkshootray(phase,delta,h,model);
end; % if nargin==6


%%% reduce p, angles and dists to those entries where p is not NaN
nonnans=find(~isnan(p));
if isempty(nonnans)
    %%% phase does not exist at this distance/depth combination
    tt=NaN;
else
    %%% phase exists, evaluate it!
    p=p(nonnans);
    angles=angles(nonnans);
    dists=dists(nonnans);

    %%% compute traveltimes for each of these ray parameters
    if ~isnan(p)
       anz=length(p); % so many solutions
       tt=p*0+NaN; % an "empty" array for the times
       for indy=1:anz
          tt(indy)=mktim4p(phase,h,p(indy),model,angles(indy));
          
          %%% control output
          %disp(['MKTTIME: ' phase ': ' num2str(h) 'km, p=',...
          %       num2str(p(indy)) ', a=' num2str(angles(indy)),...
          %       ', d=' num2str(dists(indy)),...
          %       ', derror=' num2str(dists(indy)-delta)]);
   
       end; % for indy
    else
       tt=NaN;
    end; % ~isnan
end; % if isempty(nonnans)


%%% return results
%tt=tt; already
%p=p; already