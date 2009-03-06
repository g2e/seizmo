function rayhandles=mkrayfan(phase,model,h,dangle,p5);
% mkrayfan.......create ray fan plot for a single phase
%
% call: rayhandles=mkrayfan(phase,model,h);
%       rayhandles=mkrayfan(phase,model,h,dangle);
%       rayhandles=mkrayfan(phase,model,h,colspec);
%       rayhandles=mkrayfan(phase,model,h,dangle,colspec);
%
%         phase: phase for which the ray fan is created
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
%              h: focal depth [km]
%         dangle: scalar value: take off angle resolution [deg]
%                 vector: list of take off angles [deg]
%                 DEFAULT: MKSMARTTAKEOFF will be called with resolution 1deg.
%                          usually, you may rely on this.
%         colspec: string containing a color specification.
%                  This may be either a short color name ('b', 'r', 'k', etc.),
%                  a long color name ('blue', 'red', 'black') or a RGB triple.
%                  It is important that RGB triples are given as strings ('[r g b]') also!
%                  (because otherwise it would be impossible to distinguish a RGB triple
%                  from a take off angle list)
%                  Note that this is a colorspec, not a linespec.
%
% result: rayhandles: handles to line objects of all plotted rays.
%                     Rays are drawn layer-by-layer, so a single rays consists of as many
%                     line elements as model layers are crossed. rayhandles may contain
%                     a huge number of handles...
%
% If h correpsonds to a discontinuity depth, i.e. the focus is exactly at a discontinuity,
% the above-velocity will be used for take off angles >=90deg, the below-velocity will
% be used for take off angles <90deg.
%
% The line objects are annotated:
% - the TAG field contains the phase code
% - the USERDATA field contains takeoff angle, maximum depth and epicentral
%   distance in a cell array.
%
% Martin Knapmeyer 11.06.2002, 05.07.2002, 08.09.2003, 25.09.2003,
%                  21.11.2003, 25.02.2004, 10.10.2006
%
% NOTE: This routine is based on ray parameters only. Take off angles are not considered.
%       This may result in omission of certain details!


%%% remove blanks (mey be necessary)
phase=deblank(phase);

%%% interpretation of input
nin=nargin;
switch nin
   case {3}
      %%% produce take off angle list
      dangle=1;
      angles=mksmarttakeoff(phase,model,h,dangle);
      %%% define colspec dummy
      colspec=[];
   case {4}
      if isstr(dangle)
         %%% the fourth parameter is a ColorSpec string
         colspec=str2num(dangle);
         if isempty(colspec)
            %%% colspec is empty now: dangle contains a color name.
            colspec=dangle;
         end; % if isempty(colspec)
         dangle=1;
         angles=mksmarttakeoff(phase,model,h,dangle);
      else
         %%% the fourth parameter is an angle parameter
         colspec=[];
         %%% dangle is angular resolution or list of takeoff angles
         if length(dangle)==1
            %%% angular resolution given
            angles=dangle:dangle:(180-dangle);
         else
            %%% list of angles given - make it unique MK28092006
            angles=unique(dangle(:));
         end; % if length(dangle)
      end; % if isstr(dangle)
   case {5}
      %%% p5 is a colspec
      colspec=str2num(p5);
      if isempty(colspec)
            %%% colspec is empty now: dangle contains a color name.
            colspec=p5;
         end; % if isempty(colspec)
      %%% dangle is angular resolution or list of takeoff angles
      if length(dangle)==1
         %%% angular resolution given
         angles=dangle:dangle:(180-dangle);
      else
         %%% list of angles given
         angles=dangle(:);
      end; % if length(dangle)
   otherwise
      error('MKRAYFAN: illegal number of input arguments!');
end; % swtich nin
angleanz=length(angles);

% %%% interpretation of input
% if nargin==3
%    %%% produce take off angle list
%    dangle=1;
%    angles=mksmarttakeoff(phase,model,h,dangle);
% else
%    %%% angular resolution or list of takeoff angles
%    if length(dangle)==1
%       %%% angular resolution given
%       angles=dangle:dangle:(90-dangle);
%    else
%       %%% list of angles given
%       angles=dangle(:);
%    end; % if length(dangle)
% end; % if nargin
% angleanz=length(angles);


%%% complete discontinuity list
%%% the discontinuity names defined here are referred later when searching for appropriate
%%% take off angles. Names are hard-coded in the later parts of the program! (I know it's bad style)
dz=[model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz];
dname=strvcat('Conrad','Moho','Olivine ab','Olivine bg','Olivine Perovskite','CMB','ICB',model.dname); % discontinuity names
danz=length(dz); % so many discontinuities


%%% prepare ray diagram
if ~ishold
   pos=get(gcf,'Position');
   mkraydiagram(dz,dname,model.rp,'addlabel');
   title([phase ' paths for ' model.name ', Source Depth ' num2str(h) 'km']);
   figure(gcf);
   %if pos(3)~=800
   %   set(gcf,'Position',[pos(1) 0 800 800]);
   %end; % if pos(3)~=800
end; % if ~ishold

%%% determine velocity at source depth (needed for angle-ray parm conversion)
%%% if h is a discontinuity depth, there will be two velocioty values.
%%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
focus=mkinterpmodel(model,h,'simple');
indy=find(focus.z==h);
vp=focus.vp(indy);
vs=focus.vs(indy);


%%% compute rays
rayhandles=[];
for indy=1:angleanz
    angle=angles(indy);
    %%% compute ray path
    [dist,rx,rz,rt,resp]=mkx4p(phase,h,angle,model,'angle');
    %%% add ray path to diagram - if it exists
    if (~isnan(dist))&(~isinf(dist))
       newray=mkraydiagram(rx,rz,rt,model.rp);
       rayhandles=[rayhandles newray];
       %%% apply colspec
       if ~isempty(colspec)
          set(newray,'color',colspec);
       end; % if ~isempty(colspec)
       %%% annotate in object's user data: take off angle
       set(newray,'UserData',...
           {'takeoff' num2str(angle);...
            'slo' num2str(resp);...
            'max(z)' num2str(max(rz));...
            'Delta' num2str(dist)});
    else
       %disp(['MKRAYFAN: no ray for takeoff=' num2str(angle) 'deg']);
    end; % if ~isnan & ~isinf
end; % for angle