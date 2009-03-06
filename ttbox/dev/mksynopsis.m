function syn=mksynopsis(delta,h,model,phaselist,silencer);
% mksynopsis....... create synopsis of traveltimes for a given distance and focal depth
%
% call: syn=mksynopsis(delta,h,model);
%       syn=mksynopsis(delta,h,model,phaselist);
%       syn=mksynopsis(delta,h,model,'silent');
%       syn=mksynopsis(delta,h,model,phaselist,'silent');
%
%
%           delta: epicentral distance [deg]
%               h: focal depth [km]
%           model: structure containing the velocity model
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
%          phaselist: string matrix (by strvcat) containing names of phases to evaluate
%                 DEFAULT: a comprehensive phase list will be constructed by calling
%                          mkphasecatalog('all');
%          'silent': if given, no screen output will be written (and no plot drawn)
%                    (NOTE: this string must be 'silent'!)
%
% result: syn: structure containing all computed results
%              syn.delta: the input DELTA
%              syn.h: the input H
%              syn.name: name field of the input MODEL
%              syn.names: array containing phase names
%              syn.ttimes: array containing travel times[sec]
%              syn.p: array containing ray parameters at which phases arrive
%              syn.d: array containing distances at which computed rays arrive
%              syn.a: array containing take off angles at which computet rays leave the source
%              first phases come first.
%
%              If no output variable is specified, the results will be printed to the
%              screen.
%
% Martin Knapmeyer, 04.06.2002, 05.07.2002, 28.08.2003, 26.05.2005,
%                   30.10.2006, 06.12.2006, 18.12.2006

%%% 18.12.2006 use of MKSHOOTRAY implemented

tic;

%%% prepare result
syn.delta=delta;
syn.h=h;
syn.model=model.name;
syn.names=[];
syn.ttimes=[];
syn.p=[];
syn.d=[];
syn.a=[];

%%% understand input
sliencer='off';
nin=nargin;
switch nin
   case {3}
       phaselist=mkphasecatalog('allsimple');
       silencer='off';
   case {4}
       if strcmp(lower(phaselist),'silent')==1
          phaselist=mkphasecatalog('allsimple');
          silencer='on';
       else
          silencer='off';
       end; % if strcmp(phaselist,'silencer')==1
   case {5}
       silencer='on';
   otherwise
      error('MKSYNOPSIS: illegal number of input arguments!');
end; % switch nin
phasanz=size(phaselist,1); % so many phases



%%% output format for screen
fmt='%s\t%5.2f\t%s\t%6.4f\n';



%%% does the IMODEL structure already contain the critical ray paramaters
%%% list? If not, create one! MK1812006
if ~isfield(model,'criticalrays')
    %%% imodel does not contain a cirtical rays list, create one!
    model=mkimprovemodel(model);
end; % if ~isfield(model,'criticalrays')


%%% loop through phase list
cnt=1;
for indy=1:phasanz
   phase=deblank(phaselist(indy,:));
   if strcmp(silencer,'off')==1
      disp(['MKSYNOPSIS: evaluating ' model.name ' for ' phase ' at delta=' num2str(delta)]);
   end; %if strcmp(silencer,'off')==1
   
   
   %%% use the new shooting MK18122006
   [p,a,d]=mkshootray(phase,delta,h,model);
   
   
%    if sum(isnan(p))==length(p)
%        %%% all p are NaN! retry with 360-delta distance
%        [p,a,d]=mkfindp(phase,360-delta,h,model);
%    end; % if sum(isnan(p))==length(p)
%    for indy2=1:length(p)
%      if ~isnan(p(indy2))
%          tt=mktim4p(phase,h,p(indy2),model,a(indy2));
%          syn.names=strvcat(syn.names,phase);
%          syn.ttimes(cnt)=tt;
%          syn.p(cnt)=p(indy2);
%          syn.d(cnt)=d(indy2);
%          syn.a(cnt)=a(indy2);
% 
%          cnt=cnt+1;
%      end; % if ~isnan(p(indy2)
%    end; % for indy2

       
       if sum(isnan(p))==length(p) % all p elements are NaN
          %%% use the new shooting MK18122006
          [p,a,d]=mkshootray(phase,delta,h,model);
       end; % if isempty(p)
       %if ~isnan(p(pindy))
          for indy2=1:length(p)
             if ~isnan(p(indy2))
                 tt=mktim4p(phase,h,p(indy2),model,a(indy2));
                 syn.names=strvcat(syn.names,phase);
                 syn.ttimes(cnt)=tt;
                 syn.p(cnt)=p(indy2);
                 syn.d(cnt)=d(indy2);
                 syn.a(cnt)=a(indy2);

                 cnt=cnt+1;
             end; % if ~isnan(p(indy2))
          end; % for indy2
       %end; % if ~isnan(p)
   
end; % for indy


%%% sort phases by arrival time
[syn.ttimes,sorter]=sort(syn.ttimes);
syn.names=syn.names(sorter,:);
syn.p=syn.p(sorter);
syn.a=syn.a(sorter);
syn.d=syn.d(sorter);


%%% transform distances beyond 180deg into distances below 180deg
%%% MK06122006
indies=find(syn.d>180);
syn.d(indies)=360-syn.d(indies);


%%% screen output only if no output var specified
if nargout==0
   
   if strcmp(silencer,'off')==1
      clf;
      mkraydiagram([model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz],...
                    strvcat('Conrad','Moho','Olivine ab','Olivine bg','Olivine Perovskite','CMB','ICB',model.dname),...
                    model.rp,'addlabel');
      drawnow;
      disp(['Phase arrivals for model ' model.name]);
      disp(['Delta=' num2str(delta) 'deg, z=' num2str(h) 'km']);
   end; %if strcmp(silencer,'off')==1
   for indy=1:length(syn.p)
       [dist,rx,rz,rt]=mkx4p(deblank(syn.names(indy,:)),h,syn.p(indy),model,syn.a(indy));
       syn.d(indy)=dist;
       
       if strcmp(silencer,'off')==1
          if syn.d(indy)>180
             fprintf(1,fmt,syn.names(indy,:),360-syn.d(indy),mks2hms(syn.ttimes(indy)),-syn.p(indy));
             newray=mkraydiagram(360-rx,rz,rt,model.rp);
             drawnow;
          else
             fprintf(1,fmt,syn.names(indy,:),syn.d(indy),mks2hms(syn.ttimes(indy)),syn.p(indy));
             newray=mkraydiagram(rx,rz,rt,model.rp);
             drawnow;
          end; % if syn.d
          %%% annotate in object's user data: take off angle
          set(newray,'UserData',...
              {'takeoff' num2str(syn.a(indy));...
               'slo' num2str(syn.p(indy));...
               'max(z)' num2str(max(rz));...
               'Delta' num2str(syn.d(indy));...
               'tt' num2str(syn.ttimes(indy))});
       end; %if strcmp(silencer,'off')==1
   end; % for indy
   if strcmp(silencer,'off')==1
      title([model.name ' Ray Path Synopsis for '...
             texlabel('Delta') '=' num2str(delta)  texlabel('^o')...
             ', Source Depth h=' num2str(h) 'km.']);
   end; % if strcmp(silencer,'off')==1
end; % if nargout


if strcmp(silencer,'off')==1
   toctime=toc; disp(['MKSYSNOPSIS: elapsed time ' mks2hms(toctime)]);
end; % if strcmp
