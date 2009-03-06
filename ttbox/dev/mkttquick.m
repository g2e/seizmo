function ttcurves=mkttquick(model,h,rayparm,dangle,phases,mode);
% mkttquick........quick generation of travel time curve set
%
% call: ttcurves=mkttquick(model,h,rayparm,dangle); 
%       ttcurves=mkttquick(model,h,rayparm,dangle,phases);
%       ttcurves=mkttquick(model,h,rayparm,dangle,phases,mode);
%
%         model: A structure describing the velocity distribution.
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
%         rayparm: ray parameters correpsonding to DANGLE
%                  (it is highly recommended to use output of MKSHOOTRAY as input
%                   for RAYPARM and DANGLE!)
%         dangle: if scalar: take off angle resolution [deg]
%                            takeoff angles are taken from dangle:dangle:90 [deg]
%                 if vector: list of take off angles
%                            It is recommended to take a list produced by MKSHOOTRAY here, in
%                            order to compute travel times precisely at the distances under interest!
%         phases: Phase name list. If not present or empty, a standard catalog will be used.
%                 stirng matrix like strvcat('P','SKS','ScS');
%         mode: verbosity mode. If present, no screen output will be written.
%
% result: ttcurves: structure containing traveltime curves for several phases
%                   TTCURVES consists of the following fields:
%                   .name    : model.name as given by the user
%                   .dangle  : the DANGLE goiven by the user
%                   .h       : focal depth as given by the user
%                   .anz     : number opf phases contained in structure (==length of TTC)
%                   .ttc     : sub-structure of arrays containing the data itself
%                              ttc is an array. each element of ttc describes the travel time curve
%                              of a single phase.
%                   .ttc(i).p: name of the i-th phase
%                   .ttc(i).d: epicentral distances for the i-th travel time curve [deg]
%                   .ttc(i).rayp: ray parameter for the i-th travel time curve [s/deg]
%                   .ttc(i).angle: take off angles for the i-th curve [deg]
%                   .ttc(i).t: travel time for the i-th phase
%
% based on MKTTCURVES, buit without the plotting overhead
%
% Martin Knapmeyer, 23.09.2003, 08.10.2003, 20.10.2002, 29.10.2003, 05.02.2007

%tic;

% useful
radian=pi/180;

%%% verbosity
if nargin==6
   silent=1;
else
   silent=0;
end; % if nargin

%%% prepare result structure
panz=size(phases,1); % so many phases
ttcurves.name=model.name;
ttcurves.dangle=dangle;
ttcurves.h=h;
ttcurves.anz=panz;
ttcurves.ttc=repmat(struct('p',[],'d',[],'rayp',[],'angle',[],'t',[]),ttcurves.anz,1);


%%% compute travel time curves
if ~silent
   disp(['MKTTQUICK: computing traveltime for: ']);
end; % if ~silent

%%% loop over all phases
for indy=1:panz
   cnt=1;
   tt=[];
   dist=[];
   p=[];
   rayp=[];
   pendep=[];
   takeoffangle=[];
   phase=deblank(phases(indy,:));
   if prod(size(dangle))==1
      angles=dangle:dangle:90;
   else
      angles=dangle;
   end; % if prod(size(dangle))
   
   %%% say what we do
   if ~silent
      fprintf('%s \n',phase);
   end; % if ~silent
       
   %%% loop over all take off angles    
   for anglecnt=1:length(angles)
      %%% compute only if a ray parameter exists
      if (anglecnt<=length(rayparm))
         if ~isnan(rayparm(anglecnt))
            angle=angles(anglecnt);
            switch lower(phase(1))
               case {'p'}
                    p(cnt)=rayparm(anglecnt);
               case {'s'}
                    p(cnt)=rayparm(anglecnt);
               otherwise
                    error(['MKTTQUICK: unknown phase name ' upper(phase)]);
            end; % switch lower(phase
            %disp(num2str([angle p]));
            %%% ray path end epicentral distance
            [dst,prx,prz,pty]=mkx4p(phase,h,p(cnt),model,angles(cnt));%dist(cnt)=mkx4p(phase,h,p(cnt),model);
            dist(cnt)=dst;
            rayp(cnt)=p(cnt);
            takeoffangle(cnt)=angles(cnt);
            %%% travel time
            tt(cnt)=mktim4p(phase,h,p(cnt),model,angles(cnt));
            pendep(cnt)=max(prz);
            %%% increment counter  
            cnt=cnt+1;
         else
            %%% NaN result for NaN rayparm!
            dist(cnt)=NaN; % MK24102003
            p(cnt)=NaN;
            tt(cnt)=NaN;
            pendep(cnt)=NaN;
            takeoffangle(cnt)=NaN;
            %%% increment counter  
            cnt=cnt+1;
         end; % if ~isnan(rayparm(anglecnt)
      end; % if (anglecnt<length(rayparm))
   end; % for angle
       
   %%% fill out ttcurves.ttc
   ttcurves.ttc(indy).p=phase;
   ttcurves.ttc(indy).d=dist;
   ttcurves.ttc(indy).rayp=rayp;
   ttcurves.ttc(indy).angle=takeoffangle;
   ttcurves.ttc(indy).t=tt;
       
end; % for indy

%toctime=toc; disp(['MKTTQUICK: elapsed time ' mks2hms(toctime)]);