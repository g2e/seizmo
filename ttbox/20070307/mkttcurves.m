function ttcurves=mkttcurves(model,h,dangle,phases,mode);
% mkttcurves........generate travel time curve set
%
% call: ttcurves=mkttcurves(model,h);
%       ttcurves=mkttcurves(model,h,dangle);
%       ttcurves=mkttcurves(model,h,dangle,phases);
%       ttcurves=mkttcurves(model,h,[],phases);
%       ttcurves=mkttcurves(model,h,dangle,phases,mode);
%       ttcurves=mkttcurves(model,h,[],phases,mode);
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
%         dangle: scalar value: take off angle resolution [deg]
%                               takeoff angles are taken from dangle:dangle:90 [deg]
%                 vector: list of take off angles [deg] to be used
%                         (will be used for ALL phases!)
%                 DEFAULT: if no DANGLE is given, or if DANGLE is empty,
%                          take off angles will be determined by MKSMARTTAKEOFF
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
%                   .ttc(i).d: epicentral distances for the i-th travel time curve
%                   .ttc(i).rayp: ray parameter for the i-th tarvel time curve
%                   .ttc(i).angle: take off angles for i-th travel time curve
%                   .ttc(i).t: travel time for the i-th phase
%
% based on the travel time curve plotting of MKANALYZEMODEL
%
% Martin Knapmeyer 03.07.2002, 05.07.2002, 29.08.2003, 20.10.2003, 29.10.2003
%
% BUGS: cannot handle mixed phases with P and S legs.

tic;

% useful
radian=pi/180;


%make dangle a vector
if nargin==2
   dangle=[];
else
   dangle=dangle(:);
end; % if nargin


%%% verbosity
if nargin==5
   silent=1;
else
   silent=0;
end; % if nargin


%%% control for additional plots
%%% same meaning as in MKANALYZEMODEL, but makefigs(3) is ignored here
if nargout==0
    makefigs=[0 0 1 1 1 1];
else
    makefigs=[0 0 0 0 0 0];
end; % if nargout

%%% phase catalog
%%% when adding new phases, please do not forget to add them in MKFINDENRTY, too.
if nargin<=3
   phases=[];
end; % if nargin
if (isempty(phases))
   disp('MKTTCURVES: using default phase catalog');
	sphases=[]; %mkphasecatalog('s');
	sanz=size(sphases,1); % so many s phases
	pphases=mkphasecatalog('allsimple');
	panz=size(pphases,1); % so many p phases
else
   pphases=phases; %[];
   sphases=[];
   sanz=size(sphases,1); % so many s phases
   panz=size(pphases,1); % so many p phases
end; % if nargin


%%% prepare result structure
ttcurves.name=model.name;
ttcurves.dangle=dangle;
ttcurves.h=h;
ttcurves.anz=sanz+panz;
ttcurves.ttc=repmat(struct('p',[],'d',[],'rayp',[],'angle',[],'t',[]),ttcurves.anz,1);

%%% colortable for coloring traveltime curves
stretch=2; % add this many dummy colors to anhance contrast
rgb=jet(sanz+panz+stretch);
rgb=rgb((sanz+panz+stretch):-1:1,:);

%%% prepare traveltime diagram
if makefigs(3)==1
	%disp('MKTTCURVES: prepare traveltime curve plot');
	figure(3);
	clf;
	%subplot(1,2,1);
	%mkttdiagram(60);
	title([model.name ' Travel Times, Source Depth ' num2str(h) 'km']);
	pos=get(gcf,'Position');
	%set(gcf,'Position',[pos(1) 0 650 900]);
end; % if makefigs(3)



%%% prepare ray diagram to plot nearly ANY generated ray
if makefigs(2)~=0
   disp('MKTTCURVES: prepare ray fan plot');
   dz=[model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz];
   dname=strvcat('Conrad','Moho','Olivine ab','Olivine bg','Olivine Perovskite','CMB','ICB',model.dname); % discontinuity names
   danz=length(dz); % so many discontinuities
   figure(2);
   clf;
   subplot(2,1,1);
   mkraydiagram(dz,dname,model.rp,'addlabel');
   title([model.name ', Source Depth ' num2str(h) 'km']);
   subplot(2,1,2);
   mkraydiagram(dz,dname,model.rp,'addlabel');
   pos=get(gcf,'Position');
   set(gcf,'Position',[pos(1) 0 650 900]);
   %set(gcf,'Renderer','zbuffer');
end; % if makefigs(2)

%%% prepare TT/Epidist by p diagram
if makefigs(4)~=0
   figure(4);
   clf;
end; % if makefigs(4)~=0

%%% prepare penetration depth diagram
if makefigs(5)~=0
   figure(5);
   clf;
end; % if makefigs(5)~=0

%%% prepare tau(p) diagram
if makefigs(6)~=0
   figure(6);
   clf;
end; % if makefigs(5)~=0


%%% compute velocities at focal depth
%%% needed to transform takeoff angle into ray parameter
dmy=mkinterpmodel(model,h,'simple');
vph=dmy.vp; % p velocity at focus depth
vsh=dmy.vs; % S velocity at focus depth


%%% generate TT curves
handle=[]; % to collect handles to line objects from mkttdiagram
allphases=[]; % to collect names of phases plotted into tt-curve
%%% loop through all phases
if ~silent
   disp(['MKTTCURVES: computing traveltime for: ']);
end; % if ~silent
maxtt=-inf; % to find maximum travel time in plot
if ~isempty(pphases)
	for indy=1:panz
       cnt=1;
       tt=[];
       dist=[];
       rayp=[];
       p=[];
       pendep=[];
       phase=deblank(pphases(indy,:));
       if length(dangle)==1
          if makefigs(2)==1
               angles=(dangle+indy*dangle/panz):dangle:(180-dangle); % trick: adding dangle/panz gives better ray separation
                                                            % in the optional ray fan plots.
          else
               if h==0
                  angles=dangle:dangle:(90-dangle);
               else
                  angles=dangle:dangle:(180-dangle);
               end; % if h==0
          end; % if makefigs
       else
          if isempty(dangle)
             angles=mksmarttakeoff(phase,model,h,0.5);
          else
             angles=dangle;
          end; % if isempty(dangle)
       end; % if length(dangle)
       
       if ~silent
          fprintf('%s ',phase);
       end; % if ~silent
       cspec=rgb(sanz+panz+stretch-indy,:); % color specification
       if (makefigs(2)==1)|(makefigs(2)==3)
          figure(2);
          subplot(2,1,1);
       end; % if makefigs(2)
        
       for cnt=1:length(angles)
          %%% convert take off angle into ray parameter, depending
          %%% on wave type of first leg
          switch lower(phase(1))
             case {'s'}
                vstart=vsh;
             case {'p'}
                vstart=vph;
             otherwise
                error(['MKTTCURVES: unknown phase name ' upper(phase)]);
          end; % switch
          p(cnt)=radian*sin(angles(cnt)*radian)*(model.rp-h)/vstart;
          %disp(num2str([angle p]));
          [dst,prx,prz,pty]=mkx4p(phase,h,p(cnt),model,angles(cnt));%dist(cnt)=mkx4p(phase,h,p(cnt),model);
          
          %%% store result
          dist(cnt)=dst;
          rayp(cnt)=p(cnt);
          angle(cnt)=angles(cnt);
          tt(cnt)=mktim4p(phase,h,p(cnt),model,angles(cnt));
          pendep(cnt)=max(prz);
          
          if (~isinf(tt(cnt)))&(~isnan(tt(cnt)))
             maxtt=max(maxtt,tt(cnt));
          end; % if ~isinf
          
          
          %%% plot rays - ALL rays!
           if makefigs(2)==1
              figure(2);
              subplot(2,1,1);
              rayhandle=mkraydiagram(prx,prz,pty,model.rp);
              %%% color ray according as in travel time curve plot etc
              set(rayhandle,'Color',cspec);
              %%% annotate in object's user data: take off angle
              set(rayhandle,'UserData',...
                  {'takeoff' num2str(angles(cnt));...
                    'slo' num2str(p(cnt));...
                    'max(z)' num2str(max(prz));...
                    'Delta' num2str(dst);...
                    'tt' num2str(tt(cnt))});
           end; % if makefigs 
       
       end; % for cnt
       
       %%% fill out ttcurves.ttc
       ttcurves.ttc(indy).p=phase;
       ttcurves.ttc(indy).d=dist;
       ttcurves.ttc(indy).rayp=rayp;
       ttcurves.ttc(indy).angle=angle;
       ttcurves.ttc(indy).t=tt;
          
       %%% plot
       if makefigs(3)==1
          figure(3);
          handle=[handle mkttdiagram(dist,tt,'-')];
          allphases=strvcat(allphases,phase);
          set(handle(length(handle)),'Color',cspec);
          set(handle(length(handle)),'LineWidth',0.1);
          set(handle(length(handle)),'Tag',['phase ' phase]);
       end; % if makefigs
       
       
       
      
       %%% plot t(p) and dist(p)
       %%% but before plotting, transform distance from infinite angle range to
       %%% 0...180 with reflection at 180
       dist=mod(dist,360);
       distindy=find(dist>180);
       dist(distindy)=360-dist(distindy);
       %%% and now plot it
       if makefigs(4)
          figure(4);
          subplot(2,1,1); 
          hold on
          handle2=plot(rayp,tt,'-');
          set(handle2,'Color',cspec);
          set(handle2,'LineWidth',0.1);
          set(handle2,'Tag',['phase ' phase]);
          hold off
          subplot(2,1,2);
          hold on
          handle2=plot(rayp,dist,'-');
          set(handle2,'Color',cspec);
          set(handle2,'LineWidth',0.1);
          set(handle2,'Tag',['phase ' phase]);
          hold off
       end; % if makefig(4)
       %%% plot penetration depth as function of epicentral distance
       if makefigs(5)
           figure(5);
           hold on
           pdh=plot(dist,pendep,'-');
           hold off
           set(pdh,'Color',cspec);
           set(pdh,'LineWidth',0.1);
           set(pdh,'Tag',['phase ' phase]);
           textindy=find(~isnan(dist));
           if ~isempty(textindy)
              textx=dist(textindy(1));
              texty=pendep(textindy(1));
              text(textx,texty,phase);
           end; % if ~isempty(textindy)
       end; % if makefigs(5)
       if makefigs(6)
          figure(6);
          hold on
          taup=plot(rayp,tt-rayp.*dist);
          hold off
          set(taup,'Color',cspec);
          set(taup,'LineWidth',0.1);
          set(taup,'Tag',['phase ' phase]);
       end; % if makefigs(6)
       %figure(3);
       
	end; % for panz
end; % if ~isempty pphases
if ~silent
    fprintf('\n');
end; % if ~silent




%%% generate a legend 
if makefigs(3)==1
	figure(3);
   %%%% find maximum travel time in minutes, rounded to 10min
   maxtt=ceil((maxtt/60)/5)*5;
   hold on
   mkttdiagram(maxtt);
   hold off
	legend(handle,allphases,2,'location','NorthEastOutside');
end; % if makefigs(3)

if makefigs(4)==1
   figure(4);
   subplot(2,1,1);
   xlabel('Ray Parameter [s/deg]');
   ylabel('Travel Time [sec]');
   box on
   grid on
   subplot(2,1,2);
   xlabel('Ray Parameter [s/deg]');
   ylabel('Epicentral Distance [deg]');
   box on
   grid on
end; % if makefigs(4)

if makefigs(5)==1
   figure(5);
   xlabel('Epicentral Distance [deg]');
   ylabel('Penetration Depth [km]');
   box on
   grid on
end; % if makefigs(5)

if makefigs(6)
   figure(6);
   xlabel('Ray Parameter [s/deg]');
   ylabel('\tau(p) [s]');
   box on
   grid on
end; % if makefigs(6)


toctime=toc; disp(['MKTTCURVES: elapsed time ' mks2hms(toctime)]);
