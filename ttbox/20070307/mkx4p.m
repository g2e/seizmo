function [dist,segx,segz,segtyp,resp]=mkx4p(phase,h,p,model,p5);
% mkx4p.........compute epicentral by ray parameter IN SPHERICAL EARTH
%
% call: [dist,segx,segz,segtyp,resp]=mkx4p(phase,h,p,model);
%       [dist,segx,segz,segtyp,resp]=mkx4p(phase,h,p,model,anglemode);
%       [dist,segx,segz,segtyp,resp]=mkx4p(phase,h,p,model,takeoff);
%
%       phase: string containing seismic phase name, like 'P', 'S', 'ScS', 'PKPdf',...
%              The routione does not test whether your ray partameter makes sense for
%              the phase tyou desire!
%           h: hypocentral depth [km]
%           p: ray parameter [s/deg]
%              will be transformed to the appropriate value for flat earth internally.
%              OR
%              take off angle [deg]
%              The interpretation of this parameter depends on the ANGLEMODE parameter!
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
%         anglemode: keyword to define the way input parameter P is interpreted:
%                 'angle': P is an angle and has to be translated into a ray parameter
%                 'rayparm': P is a ray parameter
%                 DEFAULT: 'rayparm'
%         takeoff: takeoff angle at source [deg]
%                  
%
%          Input parameters have to be in spherical world and are transformed 
%          to flat world internally.
%
% result: dist: epicentral distance covered by PHASE for ray parameter P, from focus to surface [deg]
%               NaN in shadow zones or whenever problems occur.
%               INF will be returned if the ray vertex is in the lowermost layer.
%               Rays that come too close to the center, where v=inf, cannot be computed
%               accurately. To avoid rays going that deep, the upper boundary of the lowermost layer
%               is set to be the deepest allowed vertex.
%         segx:  horizontal distance segments traveled in each layer [deg] IN SPHERICAL EARTH
%                together with SEGZ, this gives the ray geometry
%         segz:  vertical distance segments traveled in each layer [km] IN SPHERICAL EARTH
%                these are more or less the layer interface depths. (not for the deepest element of
%                segz, which gives the ray vertex)
%                Note: _depths_ not _radii_ !
%         segtyp: describes the wave type for easch segment of the ray:
%                 if SEGTYP(i)=='P', ray segment beginning with SEGX(i),SEGZ(i) is a p-wave
%                 if SEGTYP(i)='S', ray segment beginning with SEGX(i),SEGZ(i) is a S wave.
%                 Use this to plot S-rays and P-rays with different colors or styles.
%                 Obviously, SEGTYP is one element shorter than SEGX and SEGY.
%         resp: the ray parameter used for computation, in s/deg. This is of interest if you used
%               an angle as input.
%
% computes vertex depth of ray and solves the distance integration, using MKXP.
%
%
% DOES NOT TEST P! If you enter a ray parameter that corresponds to vertical
% or horizontal take off angle, an error will occur!
%
% Martin Knapmeyer, 30.04.2002, 02.05.2002, 06.05.2002, 07.05.2002, 05.07.2002, 25.09.2003, 13.11.2003
%                   17.11.2003, 24.03.2004, 27.08.2004, 16.02.2006,
%                   20.02.2006, 27.06.2006, 28.06.2006, 03.07.2006,
%                   04.10.2006, 11.10.2006, 12.10.2006, 19.10.2006,
%                   24.10.2006, 14.11.2006, 16.11.2006, 17.11.2006,
%                   04.12.2006
%
% LIMITATIONS: rays that go through the deepest layerare not allowed due to 
%              array indexing problems. But this does not matter, since rays close
%              to the planets center cannot be evaluated accurately.

% 13.11.2003: new calling sequence for MKRAYDEPTH, containing additional parm h
% 17.11.2003: zprecision, takeoff evaluation
% 24.02.2004: recursive handling of depth phases
% 23.03.2004: new calling sequence of MKRAYDEPTH
% 27.08.2004: SP nonsense repaired, check for complex takeoff
% 16.02.2006: SKS SEGX/SEGZ dioscrepancy repaired
% 20.02.2006: Handling of surface reflections of arbitrary degree (multiples)
% 27.06.2006: depth phase handling adapted to new version of MKRAYP2ANGLE
% 28.06.2006: generation of ScS in S code supressed
% 30.06.2006: vertex-below-CMB test blurred: vdep-model.cmb>zprecision
% 03.07.2006: angleprecision introduced.
% 04.10.2006: K leg vertext depth of SKS computed for dummy source at CMB
% 11.10.2006: takeoff==90 requires special handling in P, PP, PPP as well
%             as in S, SS, SSS, 
% 12.10.2006: takeoff==NaN caught, no SP or PS construction if second leg is NaN.
% 19.10.2006: depth phase reflection angle bug repaired
% 24.10.2006: P vertex bullshit in ScS removed, ray parameter too large
%             exception removed from ScS, because unnecessary, introduced
%             epsilon-blurring for CMB in ScS
% 14.11.2006: if takeoff==90, vdep is set to focal depth, overriding the
%             result of MKRAYDEPTH
% 16.11.2006: If no CMB or ICB exists, the respective fields in MODEL and
%             MODELSAV as well as CMBSPHER and ICBSHPER are filled with the
%             radius of the lowermost layer, just to have something in
%             there.
% 17.11.2006: But I should have done that correctly, yesterday.
% 04.12.2006  new calling sequence for MKRAYDEPTH implemented, and
%             detection of all discontinuities as needed by new MKRAYDEPTH

%tic;

%%% init results
dist=NaN;
segx=NaN;
segz=NaN;
segtyp=[]; % NaN; MK24022004

% a useful constant
radian=pi/180;

%%% vertical precision
%%% this value is used as identity threshold for comparisons of ray
%%% penetration depths and discontinuity depths, for exxample to suppress
%%% core reflections in PS ray generation
zprecision=1e-6; % [km]


%%% angular precision
%%% this value is used as identity threshold for comparisons of incidence
%%% angles at discontinuities: S waves than have incidence angles >90deg at
%%% the CMB are not S rays, but ScS or SKS - rays such steep will be
%%% filtered out in the S code. MK03072006
angleprecision=1e-6;

%%% understand input
%%% takeoff is either the true take off angle or -inf, the latter means
%%% that the takeoff angle is in [0...90], but not precisely known.
resp=[];
nin=nargin;
switch nin
   case {4}
       resp=p;
       takeoff=-inf;
   case {5}
       if isstr(p5)
          %%% input parm p might be an angle
          switch lower(p5)
             case {'rayparm'}
                 resp=p;
                 takeoff=-inf;
             case {'angle'}
                 %%% p is an angle and has to be transformed into a ray parameter
                 takeoff=p; % save angle for later use
                 p=mkangle2rayp(phase,h,takeoff,model);
                 resp=p;
             otherwise
                error(['MKX4P: unknown anglemode ' upper(p5)]);
          end; % switch anglemode
       else
          takeoff=p5;
          resp=p;
       end; % if isstr(p5)
   otherwise
      error('MKX4P: illegal number of input arguments!');
end; % switch nargin

%%% check for complex takeoff angles
%%% takeoff angle becomes complex when derived from a ray parameter
%%% that is not possible for the given wave type at the given
%%% velocity - then the requested phase does not exist.
%%% imagepsilon is used to distinguish numerical inaccuracy from
%%% really complex angles
imagepsilon=1e-6; % in units of takeoff angle
if abs(imag(takeoff))>imagepsilon
   %%% significant complex part due to nonexisting phase
   dist=NaN;
   return;
else
   %%% complex part was due to numerical inacc. and may be ignored
   takeoff=real(takeoff);
end; % if abs(imag(takeoff))


%%% check for NaN or inf ray parameters MK24.06.2005
%%% check for NaN takeoff angle MK12102006
if (isinf(p))|(isnan(p))|(isnan(takeoff))
    dist=NaN;
    return;
end; % if if (isinf(p))|(isnan(p))


%%% save input for recursive calls
%%% some phases are generated by recursive calls of MKX4P
%%% examples are PS and SP, which are constructed by separately
%%% computing P and S phases and then putting the results together
modelsav=model;
psav=p;
hsav=h;


%%% depth phase handling
%%% depth phases are handeled recursively: teh are recognized by the lower case
%%% 'p' or 's' at the first position of the phae string. This has to be computed
%%% separately, and what remains is a normal phase computation for 0km source.
%%% (this is a somewhat bad style)
firstchar=phase(1);
if strcmp(firstchar,'p')|strcmp(firstchar,'s')
   %%% depth phase recognized
   remphase=phase(2:end); % what comes after the upgoing leg
   
   
   %%% some tests
   if h==0
      %%% depth phases do not exist for surface foci
      error('MKX4P: Depth phases do not exist for surface foci!');
   end; % if h==0
   if (takeoff<90)&(~isinf(takeoff))
      %warning(['MKX4P: take off angle ' num2str(takeoff) ' invalid for depth phase!']);
      return;
   end; % if takeoff<90
   
   %%% compute P or S leg from source to surface (the depth leg)
   if isinf(takeoff)
      dlp=mkrayp2angle(upper(firstchar),h,p,model);
      dlp=dlp(1); % this is always the upgoing ray, discontinuity or not MK27062006
      anglemode='angle';
   else
      dlp=takeoff;
      anglemode='angle';
   end; % if isinf(takeoff)
   [dldist,dlsegx,dlsegz,dlsegtyp,dlresp]=mkx4p(upper(firstchar),h,dlp,model,anglemode);
   
   %%% compute remaining legs for surface focus
   if isinf(takeoff)
      %%% in this case, we have a ray parameter and the computation
      %%% of a surface-to-surface wave is straight forward.
      remp=p;
      anglemode='rayparm';
   else
      %%% in this case, we have a take off angle at the source, but we
      %%% *may* not have a ray parameter.
      %%% So we convert the angle into a ray parameter and proceed. 
      remp=mkangle2rayp(firstchar,h,takeoff,model); %%% corresponding ray parameter
      anglemode='rayparm';
   end; % if isinf(takeoff)
   [remdist,remsegx,remsegz,remsegtyp,remresp]=mkx4p(remphase,0,remp,model,anglemode);
   
   %%% fit everything together
   dist=dldist+remdist;
   segx=[dlsegx remsegx(2:end)+dldist];
   segz=[dlsegz remsegz(2:end)];
   segtyp=[dlsegtyp remsegtyp];
   resp=dlresp;
   %%% all computations complete.
   return;
end; % if strcmp(firstchar,'p')|strcmp(firstchar,'s')
%%% depth phase handling done.


%%% Handling of multiples MK20022006
%%% multiples are handled recursively: first, the "base" ray is computed
%%% (e.g. P for all multiples PP, PPP, PPPP etc.), and then multiple copies
%%% of this ray are concatenated to obtain the multiple ray.
%%% In accordance with the IASPEI nomenclature, multiples are denoted by
%%% adding the number of repetitions: P5 is PPPPP, SKS2 is SKSSKS, etc.
%%% MKX4P does currently not handle phases like PK5P (PKKKKKP).
%%% note that at this place of the code, no 'p' and 's' legs can be
%%% present, since depth phases are caught above! which means: no upgoing
%%% legs here, or: takeoff angle<=90!
[remphase,repetitions]=mkstriprepetitions(phase);
if repetitions>0
   %%% the original PHASE denotes a multiple
   if takeoff<=90
       %%% The ray is not upgoing
       
       %disp(['MKX4P: doing ' int2str(repetitions) ' repetitions of ' remphase]);

       %%% compute the base ray: from source to surface
       if isinf(takeoff)
          %%% in this case, we have a ray parameter and the computation
          %%% of a surface-to-surface wave is straight forward.
          basep=p;
          anglemode='rayparm';
       else
          %%% in this case, we have a take off angle at the source, but we
          %%% *may* not have a ray parameter.
          %%% So we convert the angle into a ray parameter and proceed. 
          basep=mkangle2rayp(remphase,h,takeoff,model); %%% corresponding ray parameter
          anglemode='rayparm';
       end; % if isinf(takeoff)
       [basedist,basesegx,basesegz,basesegtyp,baseresp]=mkx4p(remphase,h,basep,model,anglemode);


       %%% compute remaining legs: from surface to surface
       if isinf(takeoff)
          %%% in this case, we have a ray parameter and the computation
          %%% of a surface-to-surface wave is straight forward.
          remp=p;
          anglemode='rayparm';
       else
          %%% in this case, we have a take off angle at the source, but we
          %%% *may* not have a ray parameter.
          %%% So we convert the angle into a ray parameter and proceed. 
          remp=mkangle2rayp(remphase,h,takeoff,model); %%% corresponding ray parameter
          anglemode='rayparm';
       end; % if isinf(takeoff)
       [remdist,remsegx,remsegz,remsegtyp,remresp]=mkx4p(remphase,0,remp,model,anglemode);

       %%% append (repetitions-1) repetitions to the base ray
       dist=basedist; % total distance
       resp=baseresp;
       segx=basesegx;
       segz=basesegz;
       segtyp=basesegtyp;
       for indy=1:repetitions-1
           dist=dist+remdist;
           segx=[segx segx(end)+remsegx(2:end)];
           segz=[segz remsegz(2:end)];
           segtyp=[segtyp remsegtyp];
       end; % for indy

       %%% all computations complete.
       return;
   
   else
       
       %%% the ray is supposed to go upwards, but that's not allowed here
       segx=NaN;
       segy=NaN;
       segtyp=NaN;
       dist=NaN;
       return;
       
   end; % if takeoff<=90
   
end; % if repetitions>0


%%% if model is imcomplete, that is: deepest layer does not contain planets center, a dummy
%%% layer from the bootm of the model to the center is added to the model
if max(model.z)<model.rp
   anz=length(model.z); % number of layers
   model.z=[model.z' model.rp']';
   model.vp=[model.vp' model.vp(anz)]';
   model.vs=[model.vs' model.vs(anz)]';
   model.rho=[model.rho' model.rho(anz)]';
   model.qp=[model.qp' model.qp(anz)]';
   model.qs=[model.qs' model.qs(anz)]';
end; % if max(model.z)<model.rp


%%% copy model into local variables
%%% this is a historical artefact from a version which did not use the model-structure-array
%%% later versions should not use this, because it is a waste of memory!
vp=model.vp;
vs=model.vs;
rp=model.rp;
r=model.rp-model.z;

%%% check validity of P
if p==0
   warning(['MKX4P: p has to be non-zero']);
end; % if p==0
if isinf(p)
    %%% p==inf may occur when S waves are traced through liquid media...
    %%% which is of course nonsense. This seems to be the best way to catch
    %%% that. MK22022005
    dist=NaN;
end; % if isinf(p)

%%% transform ray parameter to s/rad
p=mkpdeg2rad(p);

%%% transform focal depth into a radius
h=rp-h;

%%% FLAT EARTH transform for velocity model
[vpflat,zflat]=mksfer2flat(vp,rp-r,rp);
vsflat=mksfer2flat(vs,rp-r,rp);

%%% FLAT EARTH transform for discontinuity depths
cmbspher=model.cmb; % save spherical CMB depth for use in SKS etc.
icbspher=model.icb; % save spherical ICB depth for use in SKIKS etc
alld=[model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz];
[dmy,alld]=mksfer2flat(alld,alld,model.rp);
model.conr=alld(1);
model.moho=alld(2);
model.d410=alld(3);
model.d520=alld(4);
model.d660=alld(5);
model.cmb=alld(6);
model.icb=alld(7);
model.dz=alld(8:length(alld));


%%% detect ALL discontinuities, used for MKRAYDEPTH MK04122006
disconradii=mkdetectdiscon(modelsav);


%%% if no CMB or ICB exists, the NaN values in model.cmb, model.icb,
%%% modelsav.cmb, modelsav.icb, cmbspher, and icbspher will cause trouble,
%%% since these values are often used to allow or suppress rays near the
%%% CMB or ICB - e.g. if numerical problems in vertex depth determination
%%% produces surious SKS rays. Therefore these variables are filled with
%%% the depth of the lowermost layer here. MK16112006, MK17112006
lowermostflat=zflat(end-1);
if isnan(model.cmb)
   model.cmb=lowermostflat;
   modelsav.cmb=modelsav.z(end-1);
   cmbspher=modelsav.z(end-1);
end; % if isnan(model.cmb)
if isnan(model.icb)
   model.icb=lowermostflat;
   modelsav.icb=modelsav.z(end-1);
   icbspher=modelsav.z(end-1);
end; % if isnan(model.icb)

%%% compute distances.
%%% here we have to distinguish between the different phases.

%warning on
%dbstop if warning
  
switch phase
   case {'PS'}
      %%% this phase is constructed recursively
      %%% a P-leg and an S-leg are computed for a surface focus and pasted together.
      %%% depth corection is contained in the P-leg.
      %%% Reflections at the CMB are suppressed.
      [pdist,psegx,psegz,psegtyp]=mkx4p('P',hsav,psav,modelsav,takeoff);
      if takeoff<=90
         %%% takeoff>90 would be pS!
         %%% compute new takeoff for surface reflected leg
         reflecttakeoff=mkrayp2angle('S',0,psav,modelsav); % MK27082004 takeoff=asin(psav*modelsav.vs(1)/modelsav.rp)/radian;
         %%% recursively compute converted leg
         [sdist,ssegx,ssegz,ssegtyp]=mkx4p('S',0,reflecttakeoff,modelsav,'angle'); %MK27082040 psav,modelsav,takeoff);
         if ~isnan(sdist) % Mk12102006
             if ((max(ssegz)-modelsav.cmb)<zprecision)&((max(psegz)-modelsav.cmb)<zprecision)
                % was: if (abs(max(ssegz)-modelsav.cmb)>zprecision)&(abs(max(psegz)-modelsav.cmb)>zprecision) MK28-62006
                %%% is not reflected at the CMB (it's not PScS or PcPS)
                dist=pdist+sdist;
                segx=[psegx psegx(end)+ssegx(2:end)]; %was: segx=[psegx max(psegx)+ssegx]; MK17.02.2006
                segz=[psegz ssegz(2:end)]; % was: segz=[psegz ssegz]; MK17.02.2006
                if (~isnan(pdist))&(~isnan(sdist))&(~isinf(pdist))&(~isinf(sdist))
                   segtyp=[psegtyp ssegtyp];
                end; % if ~isnan
             end; % if abs(max(ssegz)-modelsav.cmb)>zprecision
         else
             dist=NaN;
             return;
         end; % if ~isnan(sdist)
      else
         dist=NaN;
         return;
      end; % if abs(max(segz)-modelsav.cmb)<zprecision
   %%% end of PS
   
   case {'SP'}
      %%% this phase is constructed recursively
      %%% a P-leg and an S-leg are computed for a surface focus and pasted together.
      %%% depth corection is contained in the S-leg.
      %%% Reflections at CMB are suppressed
      [sdist,ssegx,ssegz,ssegtyp]=mkx4p('S',hsav,psav,modelsav,takeoff);
      if takeoff<=90
         %%% takeoff>90 would be sP
         %%% compute new takeoff for surface reflected leg
         reflecttakeoff=mkrayp2angle('P',0,psav,modelsav); % MK27082004 takeoff=asin(psav*modelsav.vp(1)/modelsav.rp)/radian;
         %%% recursively compute converted leg
         [pdist,psegx,psegz,psegtyp]=mkx4p('P',0,psav,modelsav,reflecttakeoff);
         if ~isnan(pdist) % Mk12102006
             %[pdist,psegx,psegz,psegtyp]=mkx4p('P',0,reflecttakeoff,modelsav,'angle'); %MK27082040 psav,modelsav,takeoff);
             if (abs(max(psegz)-modelsav.cmb)>zprecision)&(abs(max(ssegz)-modelsav.cmb)>zprecision)
                %%% is not reflected at the CMB (it's not ScSP or SPcP)
                dist=pdist+sdist;
                segx=[ssegx ssegx(end)+psegx(2:end)]; % was: segx=[ssegx max(ssegx)+psegx]; MK17.02.2006
                segz=[ssegz psegz(2:end)]; % was: segz=[ssegz psegz]; MK17.02.2006
                if (~isnan(pdist))&(~isnan(sdist))&(~isinf(pdist))&(~isinf(sdist))
                   segtyp=[ssegtyp psegtyp];
                end; % if ~isnan
             end; %if (abs(max(psegz)-modelsav.cmb)>zprecision)&
         else
             dist=NaN;
             return;
         end; % if ~isnan(pdist)
      else
         dist=NaN;
         return;
      end; % if abs(max(psegz)-modelsav.cmb)>zprecision
   %%% end of SP
   
   case{'P','PP','PPP','PKP','PKPPKP','PKIKP','PKIKPPKIKP'} %%% P, PKP etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% compute vertex depth
      vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
      if isnan(vdep)&(takeoff<90) % takeoff<90: mk26092003
         %%% NaN if no vertex exists
         %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Returning NaN.']);
         dist=NaN;
         return;
      else
         %%% vertex exists or ray is upgoing
         %%% FLAT EARTH transform for ray parameter
         p=mkpsfer2flat(p,rp);
         %%% FLAT EARTH transform for h (with dummy velocities)
         [dmyv,dmyz]=mksfer2flat([1 1],rp-h,rp);
         hflat=dmyz;
         %%% suppress core phase computing when vertex is outside the core! mk24052002
         %%% but only if takeoff angle<90: if takeoff>90, the ray goes up and it does
         %%% not matter if the vertex is in the core!
         if takeoff<=90 % MK11102006: takeoff<=90; was: takeoff<90, mk26092003
            %%% ray is downgoing
            %%% find shallowest possible vertex below hypocenter, forget the other ones
            indies=find(vdep<=h); % below or at hypocentral depth
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            %%% catch some exceptional vdep values
            if vdep<r(length(r)-1)
               %%% vertex is in deepest layer => abort.
               %disp('MKX4P: vertex too deep. Returning inf.');
               dist=inf;
               return;
            end;
            if isempty(vdep)
               %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Returning NaN.']);
               dist=NaN;
               return;
            end; % if isempty(vdep)
            %%% ready with catching, now continue with the ray
            
            %%% FLAT EARTH transform for vdep and h (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
            vdepflat=dmyz(1);
            hflat=dmyz(2);
            switch phase
               case {'PKP','PKPPKP'}
                  if isnan(model.cmb)
                    % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                    dist=NaN;
                    return;
                  end; % if isnan(model.cmb)
                  if model.cmb-vdepflat>zprecision % MK13102006, was: (vdepflat<=model.cmb) %|(vdepflat>model.icb)
                     %%% K-leg vertex is outside the core, this is not a core phase
                     %%% OR K-leg vertex is inside inner core - this is not PKP
                     %disp(['MKX4P: this ray parameter does not create a core phase!']);
                     dist=NaN;
                     return;
                  end; % if vdepflat
%                   %%% trial fix for the total reflection problem MK01122006
%                   if abs(vdep-(rp-modelsav.icb))<1e-6
%                      vdep=rp-modelsav.icb;
%                      [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
%                      vdepflat=dmyz(1);
%                   end; % if abs(vdep-modelsav.icb)<1e-6
%                   %%% end of trial fix MK01122006
               case {'PKIKP','PKIKPPKIKP'}
                  if isnan(model.cmb)|isnan(model.icb)
                    % disp(['MKX4P: CMB or ICB not defined - no ' phase]);
                    dist=NaN;
                    return;
                  end; % if isnan(model.icb)
                  if model.icb-vdepflat>zprecision % MK13102006, was: (vdepflat<=model.icb)|(isnan(model.icb))
                     %%% I-leg vertex is outside the inner core, this is not an inner core phase
                     %%% OR there is no inner core at all.
                     %disp('MKX4P: this ray parameter does not create an inner core phase!');
                     dist=NaN;
                     return;
                  end; % if vdepflat
               case{'P','PP','PPP'}
                  if (vdepflat-model.cmb>zprecision) % was: vdepflat>model.cmb MK30062006
                     %%% would penetrate the core, but the core is reserved for PKP and PKIKP
                     %%% this is simply to avoid confusion, a calculation of P inside the core
                     %%% would yield true paths and times!
                     dist=NaN;
                     return;
                  end; % if vdepflat
               otherwise
               % things are OK otherwise
            end; % inner switch phase
         end; % if takeoff<90 % mk26092003
         %%% compute distance from surface to real focus at HFLAT
         %disp('MKX4P: surface correction');
         [dist2,segx2,segz2]=mkxpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
         
         %%% result construction is different for upgoing and downgoing rays
         if takeoff<=90 % MK11102006was: takeoff<90 % mk26092003
            %%% compute distance for surface focus
            %disp('MKX4P: surface focus');
            [dist1,segx1,segz1]=mkxpsum(p,vpflat,zflat,0,vdepflat,0);
            %%% total epicentral distance of phase is twice the distance for surface focus
            %%% minus distance from surface to focal depth
            dist=mkxflat2sfer(2*dist1-dist2,rp);
            %%% ray geometry: several pieces have to be connected...
            %% first: transform anything into spherical world
            [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp); % for surface-to-surface
            [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp); % for focal depth correction
            segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
            segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
            %% second: build surface-to-surface ray
            indies=(length(segxs):-1:1);
            segx=[segxs 2*segxs(end)-segxs(indies)];% was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
            segz=[segz1s segz1s(indies)];
            %% third: remove ray segments from surface to focus if necessary
            if h~=rp
               indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
               segx=segx(indies);
               segz=segz(indies);
               %% fourth: shift for distance of focal depth correction
               segx=[0 segx-segx2s(length(segx2s))];
               segz=[rp-h segz];
            end; % if h~=0
            %%% generate PP and PPP if required
            %%% PP and PPP is the same as P - but two or three times.
            %%% AND the second and third branch do not need surface correction!
            switch phase
               case{'PP','PKPPKP','PKIKPPKIKP'}
                  indies=(length(segxs):-1:1);
                  segxss=[segxs 2*max(segxs)-segxs(indies)];
                  segzss=[segz1s segz1s(indies)];
                  segx=[segx segxss+dist];
                  segz=[segz segzss];
                  dist=dist+2*mkxflat2sfer(dist1,rp);
               case{'PPP'}
                  indies=(length(segxs):-1:1);
                  dist1=mkxflat2sfer(dist1,rp);
                  segxss=[segxs 2*max(segxs)-segxs(indies)]+dist;
                  segzss=[segz1s segz1s(indies)];
                  segx=[segx segxss segxss+2*dist1];
                  segz=[segz segzss segzss];
                  dist=dist+4*dist1;
               end; % switch phase (inner switch for PP and PPP)
            segtyp=char(zeros(1,length(segx)-1)+abs('P'));
         else
            %%% takeoff angle is larger then 90deg: ray goes from source upward
            %%% this happens only to P and depth phases (p-something, note
            %%% that for takeoff>90, PP becomes pP and is not handled here!)
            %%% remember: we have the depth correction ray in [dist2,segx2,segz2]
            %%% now we can transform this into spherical world
            if strcmp(phase,'P')
               [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
               segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
               dist=mkxflat2sfer(dist2,rp);
               segx=max(segx2s)-segx2s(end:-1:1);
               segz=segz2s(end:-1:1);
               segtyp=char(zeros(1,length(segx)-1)+abs('P'));
            else
               dist=NaN;
               return;
            end; % if strcmp
         end; % if takeoff<90 else % mk26092003
      end; % if isnan vdep else
      %%% end of P, PP, ... %%%
      
      %disp([' p=' num2str(p) ' z=' num2str(vdep,16) ' a=' num2str(takeoff,3) ' d=' num2str(dist)]);
      
      %%%% PKKP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case {'PKKP'}
         %%% this is a modified copy of the SKS/SKKS generation
         if takeoff<=90
            %%% to arrive at the CMB, you have to go DOWN, so takloff has to be
            %%% smaller than 90deg.
            %%% compute vertex depth for K-leg
            vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
            %disp(['MKX4P: vdep=' num2str(vdep)]);
            %%% find shallowest possible vertex below hypocenter %CMB, forget the other ones
            indies=find(vdep<=h); % was: (vdep<=model.rp-cmbspher); % below or at CMB
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            %%% catch some exceptional vdep values
            if isnan(model.cmb)
               % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
               dist=NaN;
               return;
            end; % if isnan(modelc,mb)
            if vdep<r(length(r)-1)
               %%% vertex is in deepest layer => abort.
               %disp('MKX4P: vertex too deep. Returning inf.');
               dist=inf;
               return;
            end; % if vdep<r()
            if isempty(vdep)
                %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below ' num2str(cmbspher) 'km. Returning NaN.']);
                dist=NaN;
                return;
            end; % if isempty(vdep)
            %disp(['MKX4P: vdep=' num2str(vdep)]);
            %%% ready with catching, now continue with the ray
            %%% FLAT EARTH transform for ray parameter
            p=mkpsfer2flat(p,rp);
            %%% FLAT EARTH transform for vdep and h (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
            vdepflat=dmyz(1);
            hflat=dmyz(2);
            %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' CMB flat:' num2str(model.cmb) ]);
            if (vdepflat<=model.cmb) %|(vdepflat>model.icb) %was: <
               %%% vertex is outside core: this is not PKKP.
               %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
               dist=NaN; % was: inf
               return;
            end; % if vdep
            %%% compute P-leg distance for surface focus
            %disp('MKX4P: compute PKKP surface-surface distance');
            [dist1,segx1,segz1]=mkxpsum(p,vpflat,zflat,0,model.cmb,0);
            %%% if ray does not reach CMB with this ray parameter, this is not PKKP!
            %disp([num2str(dist1) ' ' num2str(max(segz1)) ' ' num2str(model.cmb)]);
            %%% compute P-leg distance from surface to real focus at HFLAT
            %disp('MKX4P: compute PKKP focal depth correction');
            [dist2,segx2,segz2]=mkxpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
            %%% compute K-leg distance from CMB to vertex at VDEPFLAT
            %disp('MKX4P: compute PKKP K-Leg.');
            [kdist,ksegx,ksegz]=mkxpsum(p,vpflat,zflat,model.cmb,vdepflat,0);
            %ksegx=ksegx(2:length(ksegx)); % first element is redundant with end of P-leg % not any longer - commented out MK16.02.2006
                  
            %%% ray geometry: several pieces have to be connected...
            %% first: transform anything into spherical world
            [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
            [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
            [dmy,ksegzs]=mkflat2sfer(ksegz,ksegz,rp); % K-leg
            segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
            segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
            ksegxs=cumsum(mkxflat2sfer(ksegx,rp)); % K-leg
            %% second: build surface-to-surface ray
            indies=(length(segxs):-1:1);
            kindies=(length(ksegxs):-1:1);
            klegx=[ksegxs 2*max(ksegxs)-ksegxs(kindies)];
            klegz=[ksegzs ksegzs(kindies)];
            %%% the following switch could be used to construct phases with more than 2 K-legs
            switch phase
               case {'PKKP'}
                 segx=[segxs klegx+segxs(end) klegx+segxs(end)+klegx(end) 2*segxs(end)-segxs(indies)+2*klegx(end)];%was: segx=[segxs klegx+max(segxs) klegx+max(segxs)+max(klegx) 2*max(segxs)-segxs(indies)+2*max(klegx)]; MK20022006
                 segz=[segz1s klegz klegz segz1s(indies)];
                  %%% PKKP total epicentral distance contains an additional K-leg, of course
                 dist=mkxflat2sfer(2*dist1-dist2+4*kdist,rp);
            end; % inner switch: PKKP
            %error('stopped');
            %% third: remove ray segments from surface to focus if necessary
            if h~=rp
               indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
               segx=segx(indies);
               segz=segz(indies);
               %% fourth: shift for distance of focal depth correction
               segx=[0 segx-segx2s(length(segx2s))];
               segz=[rp-h segz];
            end; % if h~=0
            segtyp=char(zeros(1,length(segx)-1)+abs('P'));
            indies=find(segz>=cmbspher); % identify indices of K-legs
            indies=indies(1:(length(indies)));
            segtyp(indies)=char(zeros(1,length(indies))+abs('P')); % K-leg is a P wave!
            indies=find(diff(abs(segtyp))==3); %identify ends of K-legs. 3 is result of diff(abs('PS'))
            segtyp(indies)=char(zeros(1,length(indies))+abs('P')); % replacement before replaces one character too much
         else
            %%% takeoff >90
            dist=NaN;
         end; % if takeoff<=90
      %%%% end of PKKP
      
      %%%% PcP
      case{'PcP','PcPPcP'} %%% PcP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if takeoff<=90
            %%% you have to go DOWN to get reflected at the CMB!
            %%% compute max allowed ray parameter
            %%% uses spherical model!!!
            %%% (bear in mind: discontinuity depths in MODEL are for flat, velocities are for spher!)
            maxp=mkraydepthinv(model.rp-cmbspher,model.vp,model.rp-model.z);
            if p>maxp
               %disp(['MKX4P: Ray parameter too large! p=' num2str(p) ', maxp=' num2str(maxp)]);
               dist=inf;
               return;
            end; % if p>maxp
            %%% compute vertex depth
            vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
            %%% If there is no P vertex inside the core, return inf!
            %%% But if one exists, it has to be below the CMB.
            %%% find shallowest possible vertex below hypocenter, forget the other ones
            indies=find((vdep<=h)); % below or at hypocentral depth
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            %%% catch some exceptional vdep values
            if isnan(model.cmb)
               % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
               dist=NaN;
               return;
            end; % if isnan(modelc,mb)
            if isempty(vdep)
                 %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below r=' num2str(h) 'km. returning inf.']);
                 %dist=inf;
                 %return; % set inf and stop computation. MK 03.07.2002
                 vdep=mkflat2sfer(+1,+1,model.rp); % dummy value, just to have something to transform
            end; % if isempty(vdep)
            %%% ready with catching, now continue with the ray
            %%% FLAT EARTH transform for ray parameter
            p=mkpsfer2flat(p,rp);
            %%% FLAT EARTH transform for vdep and h (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
            vdepflat=dmyz(1);
            hflat=dmyz(2);
            if (vdepflat<model.cmb)
               %%% vertex is outside core: this is not PcP.
               %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
               dist=inf;
               return;
            end; % if vdep
            %%% compute distance for surface focus
            [dist1,segx1,segz1]=mkxpsum(p,vpflat,zflat,0,model.cmb,0);
            %%% compute distance from surface to real focus at HFLAT
            [dist2,segx2,segz2]=mkxpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
            switch phase
               case {'PcP'}
                  %%% total epicentral distance of phase is twice the distance for surface focus
                  %%% minus distance from surface to focal depth
                  dist=mkxflat2sfer(2*dist1-dist2,rp);
               case {'PcPPcP'}
                  %%% total distance is fdoru times the distance for surface focus
                  %%% minus distance from surface to focal depth
                  dist=mkxflat2sfer(4*dist1-dist2,rp);
            end; % inner switch
            %%% ray geometry: several pieces have to be connected...
            %% first: transform anything into spherical world
            [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
            [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
            segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
            segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
            %% second: build surface-to-surface ray
            switch phase
               case {'PcP'}
                  indies=(length(segxs):-1:1);
                  segx=[segxs 2*segxs(end)-segxs(indies)]; % was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006 
                  segz=[segz1s segz1s(indies)];
               case {'PcPPcP'}
                  indies=(length(segxs):-1:1);
                  segx=[segxs 2*segxs(end)-segxs(indies)]; % was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
                  segx=[segx segx+segx(end)]; % was: segx=[segx segx+max(segx)]; MK20022006
                  segz=[segz1s segz1s(indies) segz1s segz1s(indies)];
            end; % inner switch
            %% third: remove ray segments from surface to focus if necessary
            if h~=rp
               indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
               segx=segx(indies);
               segz=segz(indies);
               %% fourth: shift for distance of focal depth correction
               segx=[0 segx-segx2s(length(segx2s))];
               segz=[rp-h segz];
            end; % if h~=0
            segtyp=char(zeros(1,length(segx)-1)+abs('P'));
         else
            dist=NaN;
         end; % if takeoff<=90
      %%%% end of PcP
      
      
      %%%% PcS, ScP
      case{'PcS','ScP','PcSPcS','ScPScP'} %%% PcS, ScP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
         if takeoff<=90
            %%% you have to go DOWN to get reflected by the CMB!
             %%% compute max allowed ray parameter
             %%% uses spherical model!!!
             %%% (bear in mind: discontinuity depths in MODEL are for flat, velocities are for spher!)
             maxp=mkraydepthinv(model.rp-cmbspher,model.vp,model.rp-model.z);
             if p>maxp
                %disp(['MKX4P: Ray parameter too large! p=' num2str(p) ', maxp=' num2str(maxp)]);
                dist=inf;
                return;
             end; % if p>maxp
             %%% compute vertex depth
             vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
             %%% If there is no P vertex inside the core, return inf!
             %%% But if one exists, it has to be below the CMB.
             %%% find shallowest possible vertex below hypocenter, forget the other ones
             indies=find((vdep<=h)); % below or at hypocentral depth
             vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
             %%% catch some exceptional vdep values
             if isnan(model.cmb)
                % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                dist=NaN;
                return;
             end; % if isnan(modelc,mb)
             if isempty(vdep)
                  %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below r=' num2str(h) 'km. returning inf.']);
                  %dist=inf;
                  %return; % set inf and stop computation. MK 03.07.2002
                  vdep=mkflat2sfer(+1,+1,model.rp); % dummy value, just to have something to transform
             end; % if isempty(vdep)
             %%% ready with catching, now continue with the ray
             %%% FLAT EARTH transform for ray parameter
             p=mkpsfer2flat(p,rp);
             %%% FLAT EARTH transform for vdep and h (with dummy velocities)
             [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
             vdepflat=dmyz(1);
             hflat=dmyz(2);
             if (vdepflat<model.cmb)
                %%% vertex is outside core: this is not PcP.
                %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
                dist=inf;
                return;
             end; % if vdep
             %%% compute distance for surface focus
             [pdist1,psegx1,psegz1]=mkxpsum(p,vpflat,zflat,0,model.cmb,0); % P-leg
             [sdist1,ssegx1,ssegz1]=mkxpsum(p,vsflat,zflat,0,model.cmb,0); % S-leg
             %%% compute distance from surface to real focus at HFLAT - needed fopr depth correction
             switch phase
                case {'PcS','PcSPcS'}
                   [dist2,segx2,segz2]=mkxpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
                case {'ScP','ScPScP'}
                   [dist2,segx2,segz2]=mkxpsum(p,vsflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
             end; % inner switch
             switch phase
                case {'PcS','ScP'}
                   %%% total epicentral distance of phase is P-leg-distance plus S-leg distance
                   %%% minus distance from surface to focal depth
                   dist=mkxflat2sfer(pdist1+sdist1-dist2,rp);
                case {'PcSPcS','ScPScP'}
                   %%% total epicentral distance of phase is twice the P-leg distance plus twice th S-leg distance
                   %%% minus distance from surface to focal depth
                   dist=mkxflat2sfer(2*pdist1+2*sdist1-dist2,rp);
             end; % inner switch
             %%% ray geometry: several pieces have to be connected...
             %% first: transform anything into spherical world
             [dmy,psegz1s]=mkflat2sfer(psegz1,psegz1,rp); % P-leg
             [dmy,segz2]=mkflat2sfer(segz2,segz2,rp); % depth correction
             psegxs=cumsum(mkxflat2sfer(psegx1,rp));  % for surface-to-surface P-leg
             [dmy,ssegz1s]=mkflat2sfer(ssegz1,ssegz1,rp); % S-leg
             psegxs=cumsum(mkxflat2sfer(psegx1,rp));  % for surface-to-surface P-leg
             ssegxs=cumsum(mkxflat2sfer(ssegx1,rp));  % for surface-to-surface S-leg
             segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
             %% second: build surface-to-surface ray
             indies=(length(psegxs):-1:1);  % P ans S leg should consist of same number of pieces!
             switch phase
                case {'PcS'}
                  segx=[psegxs psegxs(end)+ssegxs(end)-ssegxs(indies(2:end))];% was: segx=[psegxs max(psegxs)+max(ssegxs)-ssegxs(indies)]; MK17022006
                  segz=[psegz1s ssegz1s(indies(2:end))];% was: segz=[psegz1s ssegz1s(indies)]; MK17022006
                  segtyp=[char(zeros(1,length(psegxs)-1)+abs('P')) char(zeros(1,length(indies)-1)+abs('S'))];
                case {'PcSPcS'}
                  segx=[psegxs psegxs(end)+ssegxs(end)-ssegxs(indies)];% was: segx=[psegxs max(psegxs)+max(ssegxs)-ssegxs(indies)]; MK17022006
                  segx=[segx segx+segx(end)];% was: segx=[segx segx+max(segx)]; MK17022006
                  segz=[psegz1s ssegz1s(indies) psegz1s ssegz1s(indies)];
                  segtyp=[char(zeros(1,length(psegxs))+abs('P'))...
                          char(zeros(1,length(ssegxs))+abs('S'))...
                          char(zeros(1,length(psegxs))+abs('P'))...
                          char(zeros(1,length(ssegxs)-1)+abs('S'))];
                case {'ScP'}
                   segx=[ssegxs ssegxs(end)+psegxs(end)-psegxs(indies(2:end))]; % was: segx=[ssegxs max(ssegxs)+max(psegxs)-psegxs(indies)]; MK17022006
                   segz=[ssegz1s psegz1s(indies(2:end))]; % was: segz=[ssegz1s psegz1s(indies)]; MK17022006
                   segtyp=[char(zeros(1,length(ssegxs)-1)+abs('S')) char(zeros(1,length(indies)-1)+abs('P'))];
                case {'ScPScP'}
                   segx=[ssegxs ssegxs(end)+psegxs(end)-psegxs(indies)];% was: segx=[ssegxs max(ssegxs)+max(psegxs)-psegxs(indies)]; MK17022006
                   segx=[segx segx+segx(end)];% was: segx=[segx segx+max(segx)]; MK17022006
                   segz=[ssegz1s psegz1s(indies) ssegz1s psegz1s(indies)];
                   segtyp=[char(zeros(1,length(ssegxs))+abs('S'))...
                           char(zeros(1,length(psegxs))+abs('P'))...
                           char(zeros(1,length(ssegxs))+abs('S'))...
                           char(zeros(1,length(psegxs)-1)+abs('P'))];
             end; % inner switch
             %% third: remove ray segments from surface to focus if necessary
             if h~=rp
                %disp(['MKX4P: depth correction: ray segment removal']);
                segzlen=length(segz);
                switch phase
                   case {'PcS','PcSPcS'}
                      indies=[find(segz(1:length(psegz1s))>(rp-h)) (length(psegz1s)+1):length(segz)]; % these remain in ray
                      removed=segzlen-length(indies); % number of removed segments
                   case {'ScP','ScPScP'}
                      indies=[find(segz(1:length(ssegz1s))>(rp-h)) (length(ssegz1s)+1):length(segz)]; % these remain in ray
                      removed=segzlen-length(indies); % number of removed segments
                end; % inner switch
	             %%% segment removal
                segx=segx(indies);
                segz=segz(indies);
                segtyp=segtyp(removed:end); %segtyp(max(1,(removed-1)):end);
                %% fourth: shift for distance of focal depth correction
                segx=[0 segx-segx2s(length(segx2s))];
                segz=[rp-h segz];
             end; % if h~=0
         else
            dist=NaN;
         end; % if takeoff<=90
      %%%% end of PcS, ScP
      
      
      %%%% PKiKP
      case{'PKiKP'} %%% PKiKP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if takeoff<=90
            %%% you have to go DOWN to reach the inner core!
             %%% compute vertex depth
             vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
             %%% existence of vertices is more or less irrelevant for PKiKP. If no vertex exists - don't care.
             %%% But if one exists, it has to be below the ICB.
             %%% find shallowest possible vertex below hypocenter, forget the other ones
             indies=find(vdep<=h); % below or at hypocentral depth
             vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
             %%% catch some exceptional vdep values
             if isnan(model.icb)|isnan(model.cmb)
                % disp(['MKX4P: CMB or ICB not defined - no ' phase]);
                dist=NaN;
                return;
             end; % if isnan(modelc,mb)
             if isempty(vdep)
                  %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Setting dummy.']);
             %    dist=NaN;
             %    return;
                  vdep=mkflat2sfer(model.icb+1,model.icb+1,model.rp); % dummy value, just to have something to transform
             end; % if isempty(vdep)
             %%% ready with catching, now continue with the ray
             %%% FLAT EARTH transform for ray parameter
             p=mkpsfer2flat(p,rp);
             %%% FLAT EARTH transform for vdep and h (with dummy velocities)
             [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
             vdepflat=dmyz(1);
             hflat=dmyz(2);
             if model.icb-vdepflat>zprecision % MK13102006was: (vdepflat<model.icb)
                %%% vertex is outside inner core: this is not PKiKP.
                %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.icb) '): this is not ' phase '. Returning inf.']);
                dist=inf;
                return;
             end; % if vdep
             %%% compute distance for surface focus
             [dist1,segx1,segz1]=mkxpsum(p,vpflat,zflat,0,model.icb,0);
             %%% compute distance from surface to real focus at HFLAT
             [dist2,segx2,segz2]=mkxpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
             %%% total epicentral distance of phase is twice the distance for surface focus
             %%% minus distance from surface to focal depth
             dist=mkxflat2sfer(2*dist1-dist2,rp);
             %%% ray geometry: several pieces have to be connected...
             %% first: transform anything into spherical world
             [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
             [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
             segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
             segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
             %% second: build surface-to-surface ray
             indies=(length(segxs):-1:1);
             segx=[segxs 2*segxs(end)-segxs(indies)];% was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
             segz=[segz1s segz1s(indies)];
             %% third: remove ray segments from surface to focus if necessary
             if h~=rp
                indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
	            segx=segx(indies);
                segz=segz(indies);
                %% fourth: shift for distance of focal depth correction
                segx=[0 segx-segx2s(length(segx2s))];
                segz=[rp-h segz];
             end; % if h~=0
             segtyp=char(zeros(1,length(segx)-1)+abs('P'));
         else
            dist=NaN;
         end; % if takeoff<=90
      %%%% end of PKiKP
      
      
   case{'S','SS','SSS'} %%% S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% compute vertex depth
      vdep=mkraydepth(p,vs,r,h,model.rp,disconradii); % vdep is a radius!
      
%       if takeoff==90
%          %%% horizontal takeoff means that the vertex is at focal depth.
%          %%% unfortunately, MKRAYDEPTH sometimes has problems to recognize
%          %%% this for numerical reasons. MK14112006
%          vdep=h;
%       end; % if takeoff==90
      
      if isnan(vdep)&(takeoff<90) % takeoff<90: mk26092003
         %%% NaN if no vertex exists
         %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Abort.']);
         dist=NaN;
         return;
      else
         %%% vertex exists or ray is upgoing.
         %%% FLAT EARTH transform for ray parameter
         p=mkpsfer2flat(p,rp);
         %%% FLAT EARTH transform for h (with dummy velocities)
         [dmyv,dmyz]=mksfer2flat([1 1],rp-h,rp);
         hflat=dmyz;
         if takeoff<=90 % MK11102006, was: takeoff<90, mk26092003
            %%% find shallowest possible vertex below hypocenter, forget the other ones
            indies=find(vdep<=h); % below hypocenter
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            %%% catch some exceptional vdep values
            if vdep<r(length(r)-1)
               %%% vertex is in deepest layer => abort.
               %disp('MKX4P: vertex too deep. Returning inf.');
               dist=inf;
               return;
            end;
            if isempty(vdep)
               %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Abort.']);
               dist=NaN;
               return;
            end; % if isempty(vdep)
            if vdep<rp-cmbspher
               %disp(['MKX4P: vertex inside core not allowed for ' phase. Abort.']);
               dist=NaN;
               return;
            end; % if vdep<rp-cmbspher
            
            %%%%% ScS cuppression code MK28062006
            %%% compute P wave vertex depth for the same ray parameter
            %%% ScS can be continued to SKS or SKIKS, this means that if a P
            %%% vertex with the same ray parameter exists within the core, the
            %%% ray parameter given in P is not an S, but an ScS/SKS phase and
            %%% has to be suppressed here!
            arriveangle=mkrayp2angle('S',cmbspher,psav,modelsav);
            arriveangle=arriveangle(1); % upper side angle!
            if (arriveangle-90>angleprecision) % was: (vdep==rp-cmbspher)&(arriveangle>90) MK29062006
                %%% S vertex is below the CMB and ray is not grazing!
                %%% This could become ScS if a P
                %%% vertext below the CMB exists.
                vdepp=mkraydepth(p,vp,r,h,model.rp,disconradii);
                indies=find(vdepp<=rp-cmbspher); 
                if ~isempty(indies)
                    %%% shallowest vertext below source is below CMB
                    %%% phase has a possible P vertex below the CMB -> not an S phase!
                    dist=NaN;
                    return;
                end; % if max(vdepp)>cmbspher
            end; % if vdep==rp-cmbspher
            %%%%% end of ScS suppression code MK28062006
            
            %%% ready with catching, now continue with the ray
            %%% FLAT EARTH transform for vdep and h (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
            vdepflat=dmyz(1);
            hflat=dmyz(2);
            if vdepflat>model.cmb %vdepflat-model.cmb>zprecision % was: vdepflat>model.cmb  MK30062006 % ">=" yields problems with S touching CMB. changed to ">" MK19052004
               %%% would penetrate the core, but the core is reserved for PKP and PKIKP
               %%% this is simply to avoid confusion, a calculation of P inside the core
               %%% would yield true paths and times!
               dist=NaN;
               return;
            end; % if vdepflat
         end; % if takeoff<90
         %%% compute distance from surface to real focus at HFLAT
         [dist2,segx2,segz2]=mkxpsum(p,vsflat,zflat,0,hflat,1);
         if takeoff<=90 % MK11102006, was: takeoff<90 mk26092003
            %%% compute distance for surface focus
            [dist1,segx1,segz1]=mkxpsum(p,vsflat,zflat,0,vdepflat,0);
            %%% total distance of phase is twice the time for surface focus
            %%% minus time from surface to focal depth
            dist=mkxflat2sfer(2*dist1-dist2,rp);
            %%% ray geometry: several pieces have to be connected...
            %% first: transform anything into spherical world
            [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
            [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
            segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
            segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
            %% second: build surface-to-surface ray
            indies=(length(segxs):-1:1);
            segx=[segxs 2*segxs(end)-segxs(indies)];% was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
            segz=[segz1s segz1s(indies)];
            %% third: remove ray segments from surface to focus if necessary
            if h~=rp
               indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
               segx=segx(indies);
               segz=segz(indies);
               %% fourth: shift for distance of focal depth correction
               segx=[0 segx-segx2s(length(segx2s))];
               segz=[rp-h segz];
            end; % if h~=0
            %%% generate SS and SSS if required
            %%% SS and SSS is the same as S - but two or three times.
            %%% AND the second and third branch are not surface corrected!
            switch phase
               case{'SS'}
                  indies=(length(segxs):-1:1);
                  segxss=[segxs 2*segxs(end)-segxs(indies)];% was: segxss=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
                  segzss=[segz1s segz1s(indies)];
                  segx=[segx segxss+dist];
                  segz=[segz segzss];
                  dist=dist+2*mkxflat2sfer(dist1,rp);
               case{'SSS'}
                  indies=(length(segxs):-1:1);
                  dist1=mkxflat2sfer(dist1,rp);
                  segxss=[segxs 2*segxs(end)-segxs(indies)]+dist; % was: segxss=[segxs 2*max(segxs)-segxs(indies)]+dist; MK20022006
                  segzss=[segz1s segz1s(indies)];
                  segx=[segx segxss segxss+2*dist1];
                  segz=[segz segzss segzss];
                  dist=dist+4*dist1;
            end; % switch phase (inner switch for PP and PPP)
            segtyp=char(zeros(1,length(segx)-1)+abs('S'));
         else
            %%% takeoff angle is larger then 90deg: ray goes from source upward
            %%% this happens only to S and depth phases (s-something, note that
            %%% for takeoff>90, SS becomes sS and is not handled here!)
            %%% An SS or SSS phase with first leg upgoing would be termed
            %%% sS or sSS and does not need to be handled here!
            %%% remember: we have the depth correction ray in [dist2,segx2,segz2]
            %%% now we can transform this into spherical world
            if strcmp(phase,'S')
               [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
               segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
               dist=mkxflat2sfer(dist2,rp);
               segx=segx2s(end)-segx2s(end:-1:1); %was: segx=max(segx2s)-segx2s(end:-1:1); MK20022006
               segz=segz2s(end:-1:1);
               segtyp=char(zeros(1,length(segx)-1)+abs('S'));
            else
               dist=NaN;
               return;
            end; % if strcmp
         end; % if takeoff<90
      end; % if isnan vdep
      %%% end of S, SS, ... %%%
      
            
      %%% ScS %%%
      case{'ScS','ScSScS'} %%% ScS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if takeoff<=90
            %%% you have to go DOWN to get reflected by the CMB!
             %%% compute max allowed ray parameter
             %%% uses spherical model!!!
             %%% (bear in mind: discontinuity depths in MODEL are for flat,
             %%% velocities are for spher!)
             %%% compute S wave vertex depth
             vdep=mkraydepth(p,vs,r,h,model.rp,disconradii); % vdep is a radius!
             indies=find(vdep<=h); % below or at hypocentral depth
             vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
             %%% for ScS, no S-vertex between focus and CMB is allowed
             %%% But if an S vertex exists, it has to be below the CMB.
             %%% find shallowest possible S vertex below hypocenter, forget the other ones
             if ~isnan(vdep)
                if vdep>rp-cmbspher+zprecision % MK24102006, was: vdep>rp-cmbspher
                   %disp(['MKX4P: S vertex above CMB not allowed for ' phase]);
                   dist=inf;
                   return;
                end; % if max(vdep)
             end; % if ~isnan
             %%% catch some exceptional vdep values
             if isnan(model.cmb)
                % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                dist=NaN;
                return;
             end; % if isnan(modelc,mb)
	          if isempty(vdep)
                  %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. returning inf.']);
                  %dist=inf;
                  %return; % set inf and stop computation MK 03.07.2002
                  vdep=mkflat2sfer(+1,+1,model.rp); % dummy value, just to have something to transform
             end; % if isempty(vdep)
             %%% ready with catching, now continue with the ray
             %%% FLAT EARTH transform for ray parameter
             p=mkpsfer2flat(p,rp);
             %%% FLAT EARTH transform for vdep and h (with dummy velocities)
             [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
             vdepflat=dmyz(1);
             hflat=dmyz(2);
% removed MK24102006
%              if vdepflat<model.cmb
%                 %%% vertex is outside core: this is not PcP.
%                 disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
%                 dist=inf;
%                 return;
%              end; % if vdep
             %%% compute distance for surface focus
             [dist1,segx1,segz1]=mkxpsum(p,vsflat,zflat,0,model.cmb,0);
             %%% compute distance from surface to real focus at HFLAT
             [dist2,segx2,segz2]=mkxpsum(p,vsflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
             switch phase
                 case {'ScS'}
                   %%% total epicentral distance of phase is twice the distance for surface focus
                   %%% minus distance from surface to focal depth
                   dist=mkxflat2sfer(2*dist1-dist2,rp);
                 case {'ScSScS'}
                   %%% total distance is four times the distance for surface focus
                   %%% minus distance from surface to focal depth
                   dist=mkxflat2sfer(4*dist1-dist2,rp);
             end; % switch
             %%% ray geometry: several pieces have to be connected...
             %% first: transform anything into spherical world
             [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
             [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
             segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
             segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
             %% second: build surface-to-surface ray
             switch phase
                case {'ScS'}
                   indies=(length(segxs):-1:1);
                   segx=[segxs 2*segxs(end)-segxs(indies)];% was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
                   segz=[segz1s segz1s(indies)];
                case {'ScSScS'}
                   indies=(length(segxs):-1:1);
                   segx=[segxs 2*segxs(end)-segxs(indies)]; % was: segx=[segxs 2*max(segxs)-segxs(indies)]; MK20022006
                   segx=[segx segx+segx(end)]; % was segx=[segx segx+max(segx)]; MK20022006
                   segz=[segz1s segz1s(indies) segz1s segz1s(indies)];
             end; % inner switch
             %% third: remove ray segments from surface to focus if necessary
             if h~=rp
                indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
	             segx=segx(indies);
                segz=segz(indies);
                %% fourth: shift for distance of focal depth correction
                segx=[0 segx-segx2s(length(segx2s))];
                segz=[rp-h segz];
             end; % if h~=0
             segtyp=char(zeros(1,length(segx)-1)+abs('S'));
         else
            dist=NaN;
         end; % if takeoff<=90
      %%% end of ScS %%%
      
      %%% SKS %%%
   case {'SKS','SKSSKS','SKKS','SKIKS'}
       if takeoff<=90
          %%% you have to go DOWN to reach the core!
             %%% compute vertex depth for S legs
             vdepsleg=mkraydepth(p,vs,r,h,model.rp,disconradii);
             %%% if no vertex below CMB exists this is not a core phase
             indies=find(((vdepsleg<=h)&(vdepsleg<=model.rp-cmbspher))|(isnan(vdepsleg)));
             %disp(num2str(vdepsleg(indies)));
             if isempty(indies)
                %disp(['MKX4P: S leg vertex above CMB (' num2str(vdepsleg) ') - this is not a core phase. Returning inf.']);
                dist=inf;
                return;
             end; % if isempty(indies)
             %%% compute vertex depth for K-leg
             vdep=mkraydepth(p,vp,r,model.rp-cmbspher,model.rp,disconradii); % MK04102006, was: vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
             %disp(['MKX4P: vdep=' num2str(vdep)]);
             %%% find vertex appropriate for phase
             switch phase
                case {'SKS','SKSSKS','SKKS'}
                   %%% find shallowest possible vertex below CMB, forget the other ones
                   if isnan(model.cmb)
                      % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                      dist=NaN;
                      return;
                   end; % if isnan(modelc,mb)
                   indies=find(vdep<=model.rp-cmbspher); % below or at CMB
                case {'SKIKS'}
                   %%% find shallowest possible vertex below ICB, forget the other ones
                   if isnan(model.cmb)|isnan(model.icb)
                      % disp(['MKX4P: CMB or ICB not defined - no ' phase]);
                      dist=NaN;
                      return;
                   end; % if isnan(modelc,mb)
                   indies=find(vdep<=model.rp-icbspher); % below or at ICB
             end; % inner switch
             vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
             %%% catch some exceptional vdep values
             if vdep<r(length(r)-1)
                %%% vertex is in deepest layer => abort.
                %disp('MKX4P: vertex too deep. Returning inf.');
                dist=inf;
                return;
             end; % if vdep<r()
             if isempty(vdep)
                 %disp(['MKX4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below CMB/ICB. Returning NaN.']);
                 dist=NaN;
                 return;
             end; % if isempty(vdep)
             %disp(['MKX4P: vdep=' num2str(vdep)]);
             %%% ready with catching, now continue with the ray
             %%% FLAT EARTH transform for ray parameter
             p=mkpsfer2flat(p,rp);
             %%% FLAT EARTH transform for vdep and h (with dummy velocities)
             [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
             vdepflat=dmyz(1);
             hflat=dmyz(2);
             switch phase
                case {'SKS','SKSSKS','SKKS'}
                     if vdepflat<model.cmb % was: (vdepflat<=model.cmb) MK28092006, even before: %|(vdepflat>model.icb)
                        %%% vertex is outside core: this is not a core phase.
                        %%% OR vertex is in inner core - this is not SKS
                        %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' not in outer core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
                        dist=NaN;
                        return;
                     end; % if vdep
                case{'SKIKS'}
                    if isnan(model.icb)
                       % disp(['MKX4P: Inner Core Bundary not defined - no ' phase]);
                       dist=NaN;
                       return;
                    end; % if isnan(modelc,mb)
                    if (vdepflat<=model.icb)
                        %%% vertex is outside core: this is not a core phase.
                        %%% OR vertex is in inner core - this is not SKS
                        %disp(['MKX4P: vertex depth (flat) ' num2str(vdepflat) ' not in inner core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
                        dist=NaN;
                        return;
                     end; % if vdep
             end; % inner switch
             %%% compute S-leg distance for surface focus
             %disp('MKX4P: compute SKS surface-surface distance');
             [dist1,segx1,segz1]=mkxpsum(p,vsflat,zflat,0,model.cmb,0);
             %%% compute S-leg distance from surface to real focus at HFLAT
             %disp('MKX4P: compute SKS focal depth correction');
             [dist2,segx2,segz2]=mkxpsum(p,vsflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
             %%% compute K-leg distance from CMB to vertex at VDEPFLAT
             %disp('MKX4P: compute SKS K-Leg.');
             [kdist,ksegx,ksegz]=mkxpsum(p,vpflat,zflat,model.cmb,vdepflat,0);
             %ksegx=ksegx(2:length(ksegx)); % first element is redundant with end of S-leg % no longer: commented out MK16022006
                      
             %%% ray geometry: several pieces have to be connected...
             %% first: transform anything into spherical world
             [dmy,segz1s]=mkflat2sfer(segz1,segz1,rp);
             [dmy,segz2s]=mkflat2sfer(segz2,segz2,rp);
             [dmy,ksegzs]=mkflat2sfer(ksegz,ksegz,rp); % K-leg
             segxs=cumsum(mkxflat2sfer(segx1,rp));  % for surface-to-surface
             segx2s=cumsum(mkxflat2sfer(segx2,rp)); % for focal depth correction
             ksegxs=cumsum(mkxflat2sfer(ksegx,rp)); % K-leg
             %% second: build surface-to-surface ray
             indies=(length(segxs):-1:1);
             kindies=(length(ksegxs):-1:1);
             klegx=[ksegxs 2*ksegxs(end)-ksegxs(kindies)];% was: klegx=[ksegxs 2*max(ksegxs)-ksegxs(kindies)]; MK20022006
             klegz=[ksegzs ksegzs(kindies)];
             %%% SKS and SKKS are very similar: SKKS simply contains an additional K leg.
             %%% SKSSKS is simply SKS doubled.
             %%% SKIKS is essentially the same as SKS.
             switch phase
                case {'SKS','SKIKS'}
                  segx=[segxs klegx+segxs(end) 2*segxs(end)-segxs(indies)+klegx(end)];% was: segx=[segxs klegx+max(segxs) 2*max(segxs)-segxs(indies)+max(klegx)]; MK20022006
                  segz=[segz1s klegz segz1s(indies)];
                  %%% SKS total epicentral distance of phase is twice the distance for surface focus
                  %%% minus distance from surface to focal depth
                  %%% plus twice the k-leg distance travelled in core
                  dist=mkxflat2sfer(2*dist1-dist2+2*kdist,rp);
                case {'SKKS'}
                  segx=[segxs klegx+segxs(end) klegx+segxs(end)+klegx(end) 2*segxs(end)-segxs(indies)+2*klegx(end)];% was: segx=[segxs klegx+max(segxs) klegx+max(segxs)+max(klegx) 2*max(segxs)-segxs(indies)+2*max(klegx)]; MK20022006
                  segz=[segz1s klegz klegz segz1s(indies)];
                  %%% SKKS total epicentral distance contains an additional K-leg, of course
                  dist=mkxflat2sfer(2*dist1-dist2+4*kdist,rp);
                case {'SKSSKS'}
                  segx=[segxs klegx+segxs(end) 2*segxs(end)-segxs(indies)+klegx(end)];% was: segx=[segxs klegx+max(segxs) 2*max(segxs)-segxs(indies)+max(klegx)]; MK20022006
                  segx=[segx segx+segx(end)];% was: segx=[segx segx+max(segx)]; MK20022006
                  segz=[segz1s klegz segz1s(indies)];
                  segz=[segz segz];
                  %%% SKS total epicentral distance of phase is twice the distance for SKS for surface-surface
                  %%% minus distance from surface to focal depth
                  %%% plus twice the k-leg distance travelled in core
                  dist=mkxflat2sfer(4*dist1-dist2+4*kdist,rp);
                  %error('stopped');
             end; % inner switch: SKS/SKKS/SKSSKS
             %error('stopped');
             %% third: remove ray segments from surface to focus if necessary
             if h~=rp
                indies=[find(segz(1:length(segz1s))>(rp-h)) (length(segz1s)+1):length(segz)]; % these remain in ray
                segx=segx(indies);
                segz=segz(indies);
                %% fourth: shift for distance of focal depth correction
                segx=[0 segx-segx2s(length(segx2s))];
                segz=[rp-h segz];
             end; % if h~=0
             segtyp=char(zeros(1,length(segx)-1)+abs('S'));
             indies=find(segz>cmbspher); % identify indices of K-legs
             indies=indies(1:(length(indies)));
             segtyp(indies)=char(zeros(1,length(indies))+abs('P')); % K-leg is a P wave!
             indies=find(diff(abs(segtyp))==3); %identify ends of K-legs. 3 is result of diff(abs('PS'))
             segtyp(indies)=char(zeros(1,length(indies))+abs('S')); % replacement before replaces one character too much
         else
            dist=NaN;
         end; % if takeoff<=90
      %%% end of SKS
   otherwise %%% unkown phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      error(['MKX4P: cannot handle phase ' phase '. Returning NaN.']);
      dist=NaN;
      return;
end; % switch phase

%%% when MKXP returns NaNs, the ray coordinate vectors are not correct.
%%% Then we reject the whole result.
if length(segx)~=length(segz)
   segx=NaN;
   segy=NaN;
   segtyp=NaN;
   dist=NaN;
end; % if length(segx)      


%%% return results
% results are already in their respective variables

%toctime=toc; disp(['MKX4P: elapsed time: ' num2str(toctime) 's.']);

