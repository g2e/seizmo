function [dist,depth,rayparm,angles]=mkdepthbydist(model,h,dangle,phase);
% mkdepthbydist.......plot penetration depth as function of epicentral distance
%
% call: [dist,depth,rayparm,angles]=mkdepthbydist(model,h,dangle,phase);
%
%       model: A structure describing the velocity distribution.
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
%              h: source depth [km]
%                 takeoff angles may be erroneous when h is identical with a discontinuity.
%         dangle: scalar value: take off angle resolution [deg]
%                               takeoff angles are taken from dangle:dangle:90 [deg]
%                 vector: list of take off angles [deg] to be used
%                         (will be used for ALL phases!)
%                 DEFAULT: if no DANGLE is given, or if DANGLE is empty,
%                          take off angles will be determined by MKSMARTTAKEOFF
%         phases: a single Phase name list. 
%
% result: dist: list of epicentral distances [deg]
%         depth: list of penetration depths [km]
%         rayparm: list of corresponding ray parameters [s/deg]
%         angles: list of take off angles corresponding to ray parameters [deg]
%
%         plot(dist,depth) will plot penetration depth as function of
%         distance, as done on the fly if the routine is called without
%         output parameter.
%
% If the model provides a list of critical ray parameters, these are used,
% but the routine does not produce such a list if none is given.
%
% Martin Knapmeyer, 09.01.2007


tic;

%%% init result
dist=[];
depth=[];






%%% compute velocities at focal depth
%%% needed to transform takeoff angle into ray parameter
dmy=mkinterpmodel(model,h,'simple');
vph=dmy.vp; % p velocity at focus depth
vsh=dmy.vs; % S velocity at focus depth

%%% phase specific take off angle list
angles=mksmarttakeoff(phase,model,h,dangle);
angleanz=length(angles); % number of agnles in list

%%% compute corresponding ray parameters
p=zeros(size(angles));
for anglecnt=1:angleanz
    p(anglecnt)=mkangle2rayp(phase,h,angles(anglecnt),model);
end; % for anglecnt


%%% list of critical ray parameters, if it exists
if isfield(model,'criticalrays')
   %%% list of critical ray parameters exists
   critp=[model.criticalrays.p; model.criticalrays.s]';
   crita=mkrayp2angle(phase,h,critp,vph,vsh,model.rp)';
   angles=[angles crita];
   p=[p critp critp]; %%% not that crita contains twice as many elements as critp!
else
   %%% list of critical ray parameters does not exist
end; % if isfield


%%% loop over angle list
angleanz=length(angles); % update number of agnles in list
dist=zeros(size(angles));
depth=zeros(size(angles));
for anglecnt=1:angleanz
    %p(anglecnt)=mkangle2rayp(phase,h,angles(anglecnt),model);
    [dst,prx,prz,pty]=mkx4p(phase,h,p(anglecnt),model,angles(anglecnt));
    dist(anglecnt)=dst;
    depth(anglecnt)=max(prz);
end; % for anglecnt


%%% sort by depth
[depth,sorter]=sort(depth);
dist=dist(sorter);
rayparm=p(sorter);
angles=angles(sorter);


%%% plot result, if no output arg is returned
if nargout==0
   clf;
   plot(dist,depth,'-');
   xlabel('Epicentral Distance [deg]');
   ylabel('Depth [km]');
   title(['Ray Turning Point / Reflection Depth, Phase: ' phase ', Velocity Model: ' model.name]);
   grid on
   axis([0 180 0 model.rp]);
   axis ij
end; % if nargout==0




toctime=toc; disp(['MKDEPTHBYDIST: elapsed time: ' num2str(toctime) 'sec.']);