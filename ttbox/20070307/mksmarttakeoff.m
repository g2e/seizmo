function [angles,minangle,maxangle]=mksmarttakeoff(phase,imodel,h,dangle);
% mksmarttakeoff.......compute phase dependently optimized takeoff angles
%
% call: [angles,minangle,maxangle]=mksmarttakeoff(phase,model,h,dangle);
%
%       phase: a seismological phase name
%       imodel: velocity model structure, as returned by MKIMPROVEMODEL
%               The model structure has to contain a criticalrays-field
%               (if it is not present, we create one on the fly, wasting
%               CPU time)
%       h: focal depth [km]
%       dangle: delta angle, take off angle resolution [deg]
%
%
% result: angles: optimized list of takeoff angles [deg]
%                 If PHASE is an unknown phase name, the full
%                 range psteepangle:dangle:90 is returned, where STEEPANGLE
%                 is the steepest takeoff angle allowe for P phases (not
%                 zero, because zero takeoff is impossible in flat earth)
%                 Empty in case of problems.
%
%         minangle: smallest possible tekaoff angle [deg]
%         maxangle: largest possible takeoff angle [deg]
%
% depending on the phase, it is not necessary to scan the whole
% 0...90deg range of takeoff angles when searching for ray paths
% or travel times. This routine produces an optimized angle list.
%
% Martin Knapmeyer, 23.09.2003, 02.09.2005, 21.02.2006, 26.06.2006,
%                   28.06.2006, 29.09.2006, 10.10.2006, 11.10.2006
%                   12.10.2006, 13.10.2006, 16.11.2006

% 21.02.2006: handling of arbitrary surface reflections (multiples)
% 26.06.2006: usage of criticalrays-field to determine CMB, ICB etc.,
%             use of mkrayp2angle to convert to angles.
% 28.06.2006: separated code for PS and SP, as well as for PcS, ScP
% 29.09.2006: dummy values for steepest allowed take off angles
% 10.10.2006: consider total reflection of S at CMB as limitation for SKS
% 11.10.2006: maxangle restricted to 90deg for PP, PPP, SS, SSS, since for
%             angles larger than 90deg, it would be pP, sS etc.
% 12.10.2006: determine maximum possible ray parameter and angle for surface
%             reflections, which is needed for PS and SP construction
% 12.10.2006: determine shallowest rays for P and S that go throught the
%             inner core - these rays limit *KIK* and *KJK* phases.
% 13.10.2006: the steepest possible takeoff angle is determined such that
%             PKIKP and SKIKS can only touch the innermost layer's top.
% 13.10.2006: depth phase handling using recursive calls.
% 16.11.2006: if no core is given, steepest angles are used instead of CMB
%             touching rays, when searching for steepest mantle phases.


%%% init result
angles=[];

%%% useful, as usually
radian=pi/180;
%angleeps=1e-3; % not needed any longer, MK03072006


%%%%%%%%%%%%%%%%%%%%%%% depth phase handling
%%% depth phases consist of two ray segments: a P or S wave which goes from
%%% the source upwards to the surface and an arbitray phase which starts at
%%% the surface and goes its way according to the ray parameter.
%%% The allowed take off angle range can therefore be determined using a
%%% recursive call, which analyzes a ray from a surface source, and a
%%% conversion of the corresponding ray parameters for the acutal focal
%%% depth afterwards. MK13102006

if (strcmp(phase(1),'s'))|(strcmp(phase(1),'p'))
   %%% PHASE designates a depth phase!
   
   %%% decompose phase name into its two parts
   depthleg=phase(1); % wave type of upgoing leg
   surfaceleg=phase(2:end); % what's reflected at the surface
   
   %%% determine angle range for the surface reflected part
   %%% here we can set a big dangle, since we're interested in minangle and
   %%% maxangle only
   [angles,minangle,maxangle]=mksmarttakeoff(surfaceleg,imodel,0,10);
   
   
   %%% ray parameters corresponding to minangle and maxangle
   %%% here we have to use the surface as source depth, again
   raypminangle=mkangle2rayp(surfaceleg,0,minangle,imodel.vp(1),imodel.vs(1),imodel.rp);
   raypmaxangle=mkangle2rayp(surfaceleg,0,maxangle,imodel.vp(1),imodel.vs(1),imodel.rp);
   
   %%% velocities at focal depth
   focus=mkinterpmodel(imodel,h,'simple');
   
   
%    %%% the ray parameter at depth h cannot be bigger than the local
%    %%% slowness r/v of the medium. To test for that, we first determine the
%    %%% maximum allowed ray parameter
%    switch depthleg
%        case {'p'}
%             v=focus.vp(1);
%        case {'s'}
%             v=focus.vs(1);
%    end; % switch depthleg
%    maxrayp=mkprad2deg(focus.rp/v);
   
   
   %%% takeoff angles at focal depth
   %%% here we convert the ray parameters into angles for a DEPHTLEG phase
   %%% upgoing from depth H. Since we want the upgoing ray, we can always
   %%% use the first of the possibly two velocities given at depth H.
   %%% take care for too large ray parameters.
%    if raypminangle<maxrayp
%        minangle=mkrayp2angle(upper(depthleg),h,raypminangle,focus.vp(1),focus.vs(1),focus.rp);
%    else
%        minangle=90;
%    end; % if raypminangle<maxrayp
%    if raypmaxangle<maxrayp
%        maxangle=mkrayp2angle(upper(depthleg),h,raypmaxangle,focus.vp(1),focus.vs(1),focus.rp);
%    else
%        maxangle=90;
%    end; % if raypmaxangle<maxrayp
   minangle=mkrayp2angle(depthleg,h,raypminangle,focus.vp(1),focus.vs(1),focus.rp);
   maxangle=mkrayp2angle(depthleg,h,raypmaxangle,focus.vp(1),focus.vs(1),focus.rp);
   
   %%% from the two solutions of minangle and maxangle, we have to take the
   %%% one which is >90deg. Since the correct solution might be NaN (in the
   %%% case if the ray parameter is tto big, the two angles might be e.g.
   %%% NaN for the one value and something <90 for the other one), we have
   %%% to chose the solution which is NOT smaller than 90deg.
   %%% Note that for a vector like [NaN 30], find(~(dmy<90)) is not the
   %%% same as find (dmy>=90) !!   MK18102006
   minangle=minangle(find(~(minangle<90)));
   maxangle=maxangle(find(~(maxangle<90))); % MK08112006
   
   %%% the resulting angles might be NaN if the ray parameters are bigger
   %%% than the medium's slowness at depth h. In that case, they're
   %%% replaced by 90deg.
   if isnan(minangle)
       minangle=90;
   end; % if isnan(minangle)
   if isnan(maxangle)
       maxangle=90;
   end; % if isnan(maxangle)
   
   
   %%% due to the reflection, minangle might be larger than maxangle now.
   %%% if so, we swap them.
   if minangle>maxangle
      swapper=minangle;
      minangle=maxangle;
      maxangle=swapper;
   end; % if minangle>maxangle
   
   
   %%% construct angle list from range and dangle
   angles=minangle:dangle:maxangle;
   if isnan(angles)
      angles=[];
   else
      if angles(end)~=maxangle
         angles=[angles maxangle];
      end; % if angles(end)
   end; % if isnan(angles)
   
   
   %%% we're done!
   %%% Now we just have to leave this routine and return the constructed
   %%% angles!
   return;
   
   
end; % 

%%%%%%%%%%%%%%%%%%%%%%% depth phase handling done.



%%% does the IMODEL structure already contain the critical ray paramaters
%%% list? If not, create one!
if ~isfield(imodel,'criticalrays')
    %%% imodel does not contain a cirtical rays list, create one!
    imodel=mkimprovemodel(imodel);
end; % if ~isfield(imodel,'criticalrays')



%%% velocities at focal depth
focus=mkinterpmodel(imodel,h,'simple');


%%% Is Vs==0 at focus?
%%% Hypocenters cannot be at depths where Vs==0, since liquids do not
%%% break. So if Vs(h)==0, there is no source, therefore no rays, therefore
%%% we quit. MK02092005
if focus.vs==0
   angles=[];
   return;
end; % if focus.vs==0


%%%%%%%%%%%%%%%%%%%%%%% max ray parm for surface reflection

%%% determine maximum possible ray parameter allowed for surface
%%% reflections: this is needed for the construction of SP and PS and the
%%% like, since it gives an upper limit for the takeoff angle
psurf=imodel.rp/imodel.vp(1);
ssurf=imodel.rp/imodel.vs(1);
%%% from radians to degrees
psurf=mkprad2deg(psurf);
ssurf=mkprad2deg(ssurf);
%%% corresponding surface take off angle
psurfangle=mkrayp2angle('P',0,psurf,imodel.vp(1),imodel.vs(1),imodel.rp);
ssurfangle=mkrayp2angle('S',0,ssurf,imodel.vp(1),imodel.vs(1),imodel.rp);

%%%%%%%%%%%%%%%%%%%%%%% max ray parm for surface reflection


%%%%%%%%%%%%%%%%%%%%%%% steepest allowed ray

%%% determine ray parameter for steepest allowed ray MK13102006
%%% since rays through the center are impossible and rays through depths
%%% close to the center are affected by large numerical errors, we restrict
%%% all rays to touching the innermost layer of the model (which may vary,
%%% depending on the model sampling)
%%% This steepest angle is no guarantee that the corresponding rays are
%%% free of numerical errors caused by the Flat Earth Transformation's zero
%%% radius singularity. It just guarantees that the ray does not hit the
%%% singularity but keeps at some distance.
%%% these angles will be used instead of CMB/ICB touch if CMB or ZICB are
%%% not present. MK16112006
psteepangle=0.1; % dummy value until real determination comes, MK29092006
ssteepangle=0.1; % dummy value until real determination comes, MK29092006

%%% ray parameters corresponding to depth of deepest layer top
[psteep,ssteep]=mkdepth2rayp(imodel,imodel.z(end-1));

%%% convert ray parameters into take off angles at the source
psteepangle=min(mkrayp2angle('P',h,psteep,focus.vp(1),focus.vs(1),imodel.rp));
ssteepangle=min(mkrayp2angle('S',h,ssteep,focus.vp(1),focus.vs(1),imodel.rp));

%%%%%%%%%%%%%%%%%%%%%%% steepest allowed ray done


%%%%%%%%%%%%%%%%%%%%%%% touch CMB

%%% determine ray parameter needed to touch CMB
%%% here we can assume that there is a criticalrays-list which alredy
%%% contains the solution of this problem.
%%% pcmb, scmb: P wave and S wave ray parameters to touch CMB
%%% pcmbangle,scmbangle: the corresponding take off angles
if ~isnan(imodel.cmb)
   %%% CMB is defined, so search for rays
   
   %%% identify critical rays of CMB discontinuity
   %%% (There might be several depth samples corresponding to CMB)
   indies=find(imodel.criticalrays.z==imodel.cmb);
   pcmb=imodel.criticalrays.p(indies); % this is probably a list
   scmb=imodel.criticalrays.s(indies); % this is probably a list
   
   %%% identify the non-NaN entries in the ray parameter sub list
   %%% The first non-NaN ray parameter should be the grazing ray
   %%% (the last one is the first ray _below_ the CMB)
   indies=find(~isnan(pcmb));
   pcmb=pcmb(indies(1));
   indies=find(~isnan(scmb)); % this selection hopefully is the same as for P
   scmb=scmb(indies(1));
   

   
   %%% if focus is at a discontinuity, FOCUS contains two velocities.
   %%% Which one we have to use depends on relative depth with respect to
   %%% CMB! (assuming that the second one is for below-side of
   %%% discontinuity) MK27062006
   if h<=imodel.cmb
       %%% source is above CMB (as usual...) use downgoing ray
       vp=focus.vp(end);
       vs=focus.vs(end);
   else
       %%% source is below CMB (extremely unlikely...) ue upgoing ray
       vp=focus.vp(1);
       vs=focus.vs(1);
   end; % if h<=imodel.cmb
   
   %%% convert ray parameters into takeoff angles
   pcmbangle=mkrayp2angle('P',h,pcmb,vp,vs,imodel.rp);
   scmbangle=mkrayp2angle('S',h,scmb,vp,vs,imodel.rp);
   
   %%% which of the two possible solutions do we have to consider?  
   if h<=imodel.cmb
       %%% source is above CMB (as usual...) use downgoing ray
       pcmbangle=pcmbangle(end); % downgoing!
       scmbangle=scmbangle(end); % downgoing!
   else
       %%% source is below CMB (extremely unlikely...) ue upgoing ray
       pcmbangle=pcmbangle(1); % upgoing!
       scmbangle=scmbangle(1); % upgoing!
   end; % if h<=imodel.cmb
   
   
else
   %%% CMB is undefined, so we cannot search for rays touching it
   %%% instead, we use the steepes possible rays MK16112006
   pcmb=psteep;
   scmb=ssteep;
   pcmbangle=psteepangle;
   scmbangle=ssteepangle;
end; % if ~isnan(imodel.cmb)
%%% now we know how to get to the CMB.

%%%%%%%%%%%%%%%%%%%%%%% touch CMB done.




%%%%%%%%%%%%%%%%%%%%%%% touch ICB

%%% determine ray parameter needed to touch ICB
if ~isnan(imodel.icb)
   %%% ICB is defined, so search for rays
   
   
   %%% identify critical rays of ICB discontinuity
   %%% (There might be several depth samples corresponding to ICB)
   indies=find(imodel.criticalrays.z==imodel.icb);
   picb=imodel.criticalrays.p(indies); % this is probably a list
   sicb=imodel.criticalrays.s(indies); % this is probably a list
   
   %%% identify the non-NaN entries in the ray parameter sub list
   %%% The first non-NaN ray parameter should be the grazing ray
   %%% (the last one is the first ray _below_ the ICB)
   indies=find(~isnan(picb));
   picb=picb(indies(1));
   indies=find(~isnan(sicb)); % this selection hopefully is the same as for P
   sicb=sicb(indies(1));
   
%    %%% identify critical rays of ICB discontinuity
%    indies=find(imodel.criticalrays.z==imodel.icb);
%   
%    
%    %%% highest non-NaN ray parameter is the one of grazing ray
%    %%% (smallest is for underside of discontinuity)
%    picb=max(imodel.criticalrays.p(indies));
%    sicb=max(imodel.criticalrays.s(indies));
   
   %%% if focus is at a discontinuity, FOCUS contains two velocities.
   %%% Which one we have to use depends on relative depth with respect to
   %%% ICB! (assuming that the second one is for below-side of
   %%% discontinuity) MK27062006
   if h<=imodel.icb
       %%% source is above ICB (as usual...), use downgoing ray
       vp=focus.vp(end);
       vs=focus.vs(end);
   else
       %%% source is below ICB (extremely unlikely...) use upgoing ray
       vp=focus.vp(1);
       vs=focus.vs(1);
   end; % if h<=imodel.icb
   
   %%% convert ray parameters into takeoff angles
   picbangle=mkrayp2angle('P',h,picb,vp,vs,imodel.rp);
   sicbangle=mkrayp2angle('S',h,sicb,vp,vs,imodel.rp);
   skicbangle=mkrayp2angle('S',h,picb,vp,vs,imodel.rp); % for Sk core phases
   
   %%% which of the two possible solutions do we have to consider?  
   if h<=imodel.cmb
       %%% source is above CMB (as usual...) use downgoing ray
       picbangle=picbangle(end); % downgoing!
       sicbangle=sicbangle(end); % downgoing!
       skicbangle=skicbangle(end);
   else
       %%% source is below CMB (extremely unlikely...) ue upgoing ray
       picbangle=picbangle(1); % upgoing!
       sicbangle=sicbangle(1); % upgoing!
       skicbangle=skicbangle(1);
   end; % if h<=imodel.cmb
   
   
   
else
   %%% ICB is undefined, so we cannot search for rays touching it
   %%% instead, we use the steepest possible rays MK16112006
   picb=psteep;
   sicb=ssteep;
   picbangle=psteepangle;
   sicbangle=ssteepangle;
   skicbangle=psteepangle;
end; % if ~isnan(imodel.icb)
%%% and now we also know how to get to the ICB.

%%%%%%%%%%%%%%%%%%%%%%% touch ICB done.


%%%%%%%%%%%%%%%%%%%%%%% total reflection at CMB
%%% determine ray parameter for steepest S and P wave that is totally reflected
%%% at the Core Mantle Boundary
%%% SKS is not limited by the ray which touches the CMB from above, but by
%%% the steepest ray which is totally reflected. The ray parameter of this
%%% ray has to be determined from Vs and Vp at/below the CMB.
%%% This evaluation is possible only if a CMB is defined.
if ~isnan(imodel.cmb)
    
    %%% identify depth samples defining the CMB
    cmbindies=find(imodel.z==imodel.cmb);
  
    
    %%% determine P velocity just below CMB
    %%% thre are two samples at CMB depth, th second one is the lower side
    vpcmb=imodel.vp(cmbindies);
    
    
    %%% determine the minimum ray parameter for total reflection
    %%% Note that this is the same for P and S, since only the P velocoty
    %%% below the CMB is important!
    %%% rays are reflected if their ray parameter is larger than this!
    cmbreflectp=(imodel.rp-imodel.cmb)/vpcmb(2);
    
    %%% from radians to degrees
    cmbreflectp=mkprad2deg(cmbreflectp);
    
    
    %%% convert this ray parameter into a takeoff angle
    pcmbreflectangle=mkrayp2angle('P',h,cmbreflectp,vp,vs,imodel.rp);
    scmbreflectangle=mkrayp2angle('S',h,cmbreflectp,vp,vs,imodel.rp);
    
    %%% which of the two possible solutions do we have to consider?  
    if h<=imodel.cmb
        %%% source is above CMB (as usual...) use downgoing ray
        pcmbreflectangle=pcmbreflectangle(end); % downgoing!
        scmbreflectangle=scmbreflectangle(end); % downgoing!
    else
        %%% source is below CMB (extremely unlikely...) ue upgoing ray
        pcmbreflectangle=pcmbreflectangle(1); % upgoing!
        scmbreflectangle=scmbreflectangle(1); % upgoing!
    end; % if h<=imodel.cmb
    
else
    %%% CMB is undefined, so we cannot search for rays reflected there!
    scmbreflectangle=NaN;
    pcmbreflectangle=NaN;
    cmbreflectp=NaN;
end; % if ~isnan(model.cmb)

%%%%%%%%%%%%%%%%%%%%%%% total reflection at CMB done.



%%%%%%%%%%%%%%%%%%%%%%% shallowest ray takeoff
%%% determine ray parameter for shallowest P and S rays that can go throught
%%% the inner core: These limit the inner core phases *KIK* and *KJK*
%%% towards larger ray parameters/take off angles. MK12102006
if ~isnan(imodel.icb)
    
    %%% identify depth samples that define the ICB
    icbindies=find(imodel.z==imodel.icb);
    
    %%% determine P and S velocities just below the ICB
    %%% there are two samples, the second one is the lower side of the
    %%% discontinuity
    vpicb=imodel.vp(icbindies(2));
    vsicb=imodel.vs(icbindies(2));
    
    %%% determine the ray parametes for having a vertex there
    picbshallowp=(imodel.rp-imodel.icb)/vpicb;
    sicbshallowp=(imodel.rp-imodel.icb)/vsicb;
    
    %%% from radians to degrees
    picbshallowp=mkprad2deg(picbshallowp);
    sicbshallowp=mkprad2deg(sicbshallowp);
    
    %%% convert ray parameters into takeoff angles
    picbshallowangle=mkrayp2angle('P',h,picbshallowp,vp,vs,imodel.rp);
    sicbshallowangle=mkrayp2angle('S',h,sicbshallowp,vp,vs,imodel.rp);
    skicbshallowangle=mkrayp2angle('S',h,picbshallowp,vp,vs,imodel.rp); % for K legs in S core phases
    
    %%% which of the two possible solutions do we have to consider?
    if h<=imodel.icb
        %%% source is above ICB (as usual...) use downgoing ray
        picbshallowangle=picbshallowangle(end); % downgoing!
        sicbshallowangle=sicbshallowangle(end); % downgoing!
        skicbshallowangle=skicbshallowangle(end);
    else
        %%% source is below ICB (extremely unlikely...) ue upgoing ray
        picbshallowangle=picbshallowangle(1); % upgoing!
        sicbshallowangle=sicbshallowangle(1); % upgoing!
        skicbshallowangle=skicbshallowangle(1);
    end; % if h<=imodel.cmb
    
    
else
    %%% Inner Core Boundary is not defined - ray construction impossible
    sicbshallowangle=NaN;
    picbshallowangle=NaN;
    skicbshallowangle=NaN;
    skicbshallowp=NaN;
    picbshallowp=NaN;
    sicbshallowp=NaN;
end; % if ~isnan(imodel.icb)

%%%%%%%%%%%%%%%%%%%%%%% shallowest ray take off done.






%%% surface reflections handling MK21.02.2006
%%% here it is sufficient just to strip the repeat-factor off the phase name
[phase,repetitions]=mkstriprepetitions(phase);


%%% define phase dependent angle range
switch phase
   case {'P'}        
        minangle=pcmbangle;
        if isnan(minangle)
           minangle=psteepangle;
        end; % if isnan(minangle)
        if h>0
          maxangle=180;
        else
          maxangle=90;
        end; % if h>0
   case {'PP','PPP'}        
        minangle=pcmbangle;
        if isnan(minangle)
           minangle=psteepangle;
        end; % if isnan(minangle)
        maxangle=90;  % takeoff>90 would be called pP!   
   case {'PcP','PcPPcP'}        
        minangle=psteepangle;
        maxangle=pcmbangle;
   case {'PKP','PKKP','PKPPKP'}        
        minangle=psteepangle;
        maxangle=min(pcmbreflectangle,pcmbangle);
  case {'PKiKP'}        
        minangle=psteepangle;
        maxangle=picbangle;
   case {'PKIKP','PKIKPPKIKP'}        
        minangle=psteepangle;
        maxangle=picbshallowangle;
   case {'S'}        
        minangle=scmbangle;
        if isnan(minangle)
           minangle=ssteepangle;
        end; % if isnan(minangle)
        if h>0
          maxangle=180;
        else
          maxangle=90;
        end; % if h>0
   case {'SS','SSS'}        
        minangle=scmbangle;
        if isnan(minangle)
           %%% this means: there is no core defined
           minangle=ssteepangle;
        end; % if isnan(minangle)
        maxangle=90; % takeoff>90 would be called sS!
   case {'ScS','ScSScS'}        
        minangle=ssteepangle;
        if isnan(imodel.icb)
           %%% no inner core defined, which means that ssteepangle is zero
           %%% (steepest ScS ray might be zero, but MKX4P blocks vertical rays,
           %%%  so we use something else here)
           minangle=psteepangle;
        else
           minangle=ssteepangle;
        end; % if isnan(icmodel.icb)
        maxangle=scmbangle;
   case {'SKS','SKKS','SKSSKS'}      
        if imodel.vs(end-1)==0
            %%% S waves not possible in the lowermost layers, so the
            %%% steepest possible ray is determined by P take off angles
            %%% MK08112006
            minangle=psteepangle;
        else
            %%% S waves possible, use S take off angles.
            minangle=ssteepangle;
        end; % if imodel.vs(end-1)==0
        maxangle=min(scmbreflectangle,scmbangle);
   case {'SKiKS'}        
        minangle=ssteepangle;
        maxangle=skicbangle; % it's the K leg which is important!     
   case {'SKIKS'}        
        minangle=ssteepangle;
        maxangle=skicbshallowangle; % it's the I leg which is important!
   case {'PS'}
        if ~isnan(imodel.cmb)
            %%% model has a core, steepest S ray determined by CMB
            minangle=mkrayp2angle('P',h,scmb,vp,vs,imodel.rp); %MK12102006 was: scmbangle;
        else
            %%% model does not have a core, steepest S ray determined by
            %%% lowermost layer MK16112006
            minangle=mkrayp2angle('P',h,ssteep,focus.vp(end),focus.vs(end),imodel.rp);
        end; % if ~isnan(imodel.cmb)
        minangle=minangle(end);
        if isnan(minangle)
           minangle=psteepangle;
        end; % if isnan(minangle)
        maxangle=min([psurfangle 90]); % takeoff>90 would be called pS!
   case {'SP'}
        minangle=scmbangle;
        if isnan(minangle)
           minangle=ssteepangle;
        end; % if isnan(minangle)
        maxangle=min([psurfangle 90]); % takeoff>90 would be called sP!
   case {'PcS'}
        minangle=psteepangle; % MK29092006
        maxangle=pcmbangle;
   case {'ScP'}
        minangle=psteepangle; % MK29092006
        maxangle=scmbangle;
   otherwise
       disp(['MKSMARTTAKEOFF: unknown phase ' phase ', using default.']);
       minangle=psteepangle;
       if h>0
          maxangle=180;
       else
          maxangle=90;
       end; % if h>0
end; % switch phase



%%% construct angle list from range and dangle
angles=minangle:dangle:maxangle;

if isnan(angles)
   angles=[];
else
   if angles(end)~=maxangle
      angles=[angles maxangle];
   end; % if angles(end)
end; % if isnan(angles)