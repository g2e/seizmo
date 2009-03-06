function [tt,vdep,resp]=mktim4p(phase,h,p,model,p5);
% mktim4p.........compute traveltime by ray parameter IN SPHERICAL EARTH
%
% call: [tt,vdep,resp]=mktim4p(phase,h,p,model);
%       [tt,vdep,resp]=mktim4p(phase,h,p,model,anglemode);
%       [tt,vdep,resp]=mktim4p(phase,h,p,model,takeoff);
%
%       phase: string containing seismic phase name, like 'P', 'S', 'ScS', 'PKPdf',...
%           h: hypocentral depth [km]
%           p: ray parameter [s/deg]
%              will be transformed to the appropriate value for flat earth internally.
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
%
%         anglemode: keyword to define the way input parameter P is interpreted:
%                 'angle': P is an angle and has to be translated into a ray parameter
%                 'rayparm': P is a ray parameter
%                 DEFAULT: 'rayparm'
%         takeoff: takeoff angle at source [deg]
%
%          Input parameters have to be in spherical world and are transformed 
%          to flat world internally.
%
% result:   tt: traveltime of PHASE for ray parameter P, from focus to surface [s]
%               NaN in shadow zones or whenever problems occur.
%         vdep: radius at which the vertex is [km]
%               (the deepest vertex for mixed phases like PS)
%         resp: the ray parameter used for computation. This is of interest if you used
%            an angle as input.
%
% computes vertex depth of ray and solves the traveltime integration, using MKTP.
%
% Martin Knapmeyer, 24.04.2002, 05.07.2002, 25.09.2003, 13.11.2003,
%                   17.11.2003, 24.02.2004, 27.08.2004, 21.02.2006
%                   27.06.2006, 28.06.2006, 30.06.2006, 03.07.2006,
%                   04.10.2006, 11.10.2006, 12.10.2006, 19.10.2006,
%                   24.10.2006, 16.11.2006, 17.11.2006, 04.12.2006

% 17.11.2003: zprecision and evaluation of takeoff
% 24.02.2004: recursive handling of depth phases
% 24.03.2004: new calling sequence of MKRAYDEPTH
% 27.08.2004: SP nonsense repaired, check for complex takeoff
% 21.02.2006: surface reflections (multiples) of arbitrary degree
% 27.06.2006: depth phase handling adapted to new version of MKRAYP2ANGLE
% 28.06.2006: generation of ScS in S case suppressed
% 30.06.2006: vertex-below-CMB test blurred: vdep-(rp-cmbsher)>zprecision
% 03.07.2006: angleprecision introduced.
% 04.10.2006: K leg vertext depth for SKS computed for dummy source at CMB
% 11.10.2006: takeoff==90 requires special handling in P, PP, PPP as well
%             as in S, SS, SSS, 
% 12.10.2006: takoff==NaN caught, no PS or SP construction if second leg
%             returns NaN distance
% 19.10.2006: depth phase reflection angle bug repaired
% 24.10.2006: removed P vertext bullshit from ScS, ray parameter too large
%             exception removed from ScS, because unnecessary
% 16.11.2006: If no CMB or ICB exists, the respective fields in MODEL and
%             MODELSAV as well as CMBSPHER and ICBSHPER are filled with the
%             radius of the lowermost layer, just to have something in
%             there.
% 17.11.2006: But I should have done this correctly, yesterday.
% 04.12.2006  new calling sequence for MKRAYDEPTH implemented, and
%             detection of all discontinuities as needed by new MKRAYDEPTH



%tic;


%%% init result
tt=NaN;
vdep=NaN;
resp=NaN;


%%% a useful constant
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
                error(['MKTIM4P: unknown anglemode ' upper(anglemode)]);
          end; % switch anglemode
       else
          takeoff=p5;
          resp=p;
       end; % if isstr(p5)
   otherwise
      error('MKTIM4P: illegal number of input arguments!');
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
   tt=NaN;
   return;
else
   %%% complex part was due to numerical inacc. and may be ignored
   takeoff=real(takeoff);
end; % if abs(imag(takeoff))

%%% check for NaN or inf ray parameters MK24.06.2005
%%% check for NaN takeoff angle MK12102006
if (isinf(p))|(isnan(p))|(isnan(takeoff))
    dist=NaN;
    return
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
      error('MKTIM4P: Depth phases do not exist for surface foci!');
   end; % if h==0
   if (takeoff<90)&(~isinf(takeoff))
      %warning(['MKTIM4P: take off angle ' num2str(takeoff) ' invalid for depth phase!']);
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
   [dltt,dlvdep,dlresp]=mktim4p(upper(firstchar),h,dlp,model,anglemode);
   
%    %%% compute remaining legs for surface focus
%    if isinf(takeoff)
%       remp=p;
%       anglemode='rayparm';
%    else
%       remp=180-takeoff;
%       anglemode='angle';
%    end; % if isinf(takeoff)
%    [remtt,remvdep,remresp]=mktim4p(remphase,0,remp,model,anglemode);

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
   [remtt,remvdep,remresp]=mktim4p(remphase,0,remp,model,anglemode);
   
   
   %%% fit everything together
   tt=dltt+remtt;
   vdep=max(dlvdep,remvdep);
   resp=dlresp;
   
   %%% all computations complete.
   return;
end; % 
%%% depth phase handling done.


%%% Handling of multiples MK21.02.2006
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
       %disp(['MKTIM4P: doing ' int2str(repetitions) ' repetitions of ' remphase]);

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
          basep=mkangle2rayp(firstchar,h,takeoff,model); %%% corresponding ray parameter
          anglemode='rayparm';
       end; % if isinf(takeoff)
       [basett,basevdep,baseresp]=mktim4p(remphase,h,basep,model,anglemode);


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
          remp=mkangle2rayp(firstchar,h,takeoff,model); %%% corresponding ray parameter
          anglemode='rayparm';
       end; % if isinf(takeoff)
       [remtt,remvdep,remresp]=mktim4p(remphase,0,remp,model,anglemode);

       %%% append (repetitions-1) repetitions to the base ray
       tt=basett;
       resp=baseresp;
       vdep=min(basevdep,remvdep); % it's a radius, therefore the smaller is deeper!
       for indy=1:repetitions-1
           tt=tt+remtt;
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


%%% transform ray parameter to s/rad
p=mkpdeg2rad(p);

%%% copy model into local variables
%%% this is a historical artefact from a version which did not use the model-structure-array
%%% later versions should not use this, because it is a waste of memory!
%   vp: P wave velocities at upper and lower boundary of layer [km/s]
%       VP(1): upper boundary, VP(2): lower boundary
%   vs: S wave velocities at upper and lower boundary of layer [km/s]
%       VS(1): upper boundary, VS(2): lower boundary
%    r: radii of upper and lower boundary of layer [km]
%       R(1): upper boundary, R(2): lower boundary
%   rp: planetary radius [km]
vp=model.vp;
vs=model.vs;
rp=model.rp;
r=model.rp-model.z;


%%% transform focal depth into a radius
h=rp-h;


%%% FLAT EARTH transform for velocity model
[vpflat,zflat]=mksfer2flat(vp,rp-r,rp);
vsflat=mksfer2flat(vs,rp-r,rp);


%%% FLAT EARTH transform for discontinuity depths
cmbspher=model.cmb; % for later use with SKS etc.
icbspher=model.icb; % for later use with SKIKS etc.
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
%%% the radius of the lowermost layer here. MK16112006, 17112006
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




%%% compute times.
%%% here we have to distinguish between the different phases.
  
switch phase
   case {'PS'}
      %%% this phase is constructed by recursive calls for P-leg and S-leg
      %%% results of these calls are then added
      %%% surface correction is contained in the p-leg
      %%% reflections at CMB are suppressed
      [ptt,pvdep]=mktim4p('P',hsav,psav,modelsav,takeoff);
      if takeoff<=90
         %%% takeoff>90 would be pS
         %%% compute new takeoff for surface reflected leg
         reflecttakeoff=mkrayp2angle('S',0,psav,modelsav); % MK27082004 takeoff=asin(psav*modelsav.vs(1)/modelsav.rp)/radian;
         %%% recursively compute converted leg
         [stt,svdep]=mktim4p('S',0,reflecttakeoff,modelsav,'angle'); % MK27082004 psav,modelsav,takeoff);
         if ~isnan(stt) % MK12102006
             if (abs(svdep-modelsav.cmb)>zprecision)&(abs(pvdep-modelsav.cmb)>zprecision)
                %%% ray is not reflected at CMB
                tt=ptt+stt;
                vdep=min(svdep,pvdep);
             end; % if abs(svdep-modelsav.cmb)>zprecision
         else
             tt=NaN;
         end; % if ~isnan(stt)
      else
         tt=NaN;
      end; % 
   %%% end of PS
   case {'SP'}
      %%% this phase is constructed by recursive calls for P-leg and S-leg
      %%% results of these calls are then added
      %%% surface correction is contained in the p-leg
      [stt,svdep]=mktim4p('S',hsav,psav,modelsav,takeoff);
      if takeoff<=90
         %%% takeoff>90 would be sP
         %%% compute new takeoff for surface reflected leg
         reflecttakeoff=mkrayp2angle('P',0,psav,modelsav); % MK27082004 takeoff=asin(psav*modelsav.vp(1)/modelsav.rp)/radian;
         %%% recursively compute converted leg
         [ptt,pvdep]=mktim4p('P',0,reflecttakeoff,modelsav,'angle'); %MK27082004 psav,modelsav,takeoff);
         if ~isnan(ptt) % MK12102006
             if (abs(svdep-modelsav.cmb)>zprecision)&(abs(pvdep-modelsav.cmb)>zprecision)
                tt=ptt+stt;
                vdep=min(svdep,pvdep);
             end; % if (abs(svdep-modelsav.cmb)>zprecision)&
         else
             tt=NaN;
         end; % if ~isnan(ptt)
      else
         tt=NaN;
      end; % if (abs(svdep-modelsav.cmb)>zprecision)&
   %%% end of PS
   case{'P','PP','PPP','PKP','PKPPKP','PKIKP','PKIKPPKIKP'} %%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% compute vertex depth
      vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
      if isnan(vdep)&(takeoff<90) % takeoff<90: mk26092003
         %%% NaN if no vertex exists
         %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Returning NaN.']);
         tt=NaN;
         return;
      else
         %%% vertex exists or ray is upgoing.
         %%% FLAT EARTH transform for ray parameter
         p=mkpsfer2flat(p,rp);
         %%% FLAT EARTH transform for h (with dummy velocities)
         [dmyv,dmyz]=mksfer2flat([1 1],rp-h,rp);
         hflat=dmyz;
         if takeoff<=90 % MK11102006, was: takeoff<90 % if takeoff<90 MK26092003
            %%% find shallowest possible vertex below hypocenter, forget the other ones
            indies=find(vdep<=h); % below or at hypocentral depth
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            if vdep<r(length(r)-1)
               %%% vertex is in deepest layer => abort.
               %disp('MKX4P: vertex too deep. Abort.');
               tt=inf;
               return;
            end;
            if isempty(vdep)
               %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Returning NaN.']);
               tt=NaN;
               return;
            end; % if isempty(vdep)
            %%% catch some exceptions
            switch phase
               case {'PKP','PKPPKP'}
                  if isnan(model.cmb)
                       % disp(['MKTIM4P: Core Mantle Boundary not defined - no ' phase]);
                       tt=NaN;
                       return;
                  end; % if isnan(model.cmb)
                  if vdep>=rp-cmbspher
                     %%% K-leg vertex is outside the core, this is not a core phase
                     %%% OR K-leg vertex is inside inner core - this is not PKP
                     %disp(['MKX4P: this ray parameter does not create a core phase!']);
                     tt=NaN;
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
                  if isnan(model.cmb)
                    % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                    tt=NaN;
                    return;
                  end; % if isnan(model.cmb)
                  if isnan(model.icb)
                    % disp(['MKX4P: Inner core Boundary not defined - no ' phase]);
                    tt=NaN;
                    return;
                  end; % if isnan(model.cmb)
                  if vdep>=rp-icbspher %(vdepflat<=model.icb)|(isnan(model.icb))
                     %%% I-leg vertex is outside the inner core, this is not an inner core phase
                     %%% OR there is no inner core at all.
                     %disp('MKX4P: this ray parameter does not create an inner core phase!');
                     tt=NaN;
                     return;
                  end; % if vdepflat
               case{'P','PP','PPP'}
                  if vdep-(rp-cmbspher)<zprecision % was: vdep<rp-cmbspher MK30062006
                     %%% would penetrate the core, but the core is reserved for PKP and PKIKP
                     %%% this is simply to avoid confusion, a calculation of P inside the core
                     %%% would yield true paths and times!
                     tt=NaN;
                     return;
                  end; % if vdepflat
               otherwise
                  % otherwise, things are OK
            end; % inner switch 
         end; % if takeoff<90
         %%% compute traveltime from surface to real focus at HFLAT
         tt2=mktpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
         if takeoff<=90 % MK11102006, was: takeoff<90
            %%% FLAT EARTH transform for vdep (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-vdep,rp);
            vdepflat=dmyz;
            %%% compute traveltime for surface focus
            tt1=mktpsum(p,vpflat,zflat,0,vdepflat,0);
         
            %%% total traveltime of phase is twice the time for surface focus
            %%% minus time from surface to focal depth
            tt=2*tt1-tt2;
            %%% generate PP and PPP if required
            %%% PP and PPP is the same as P - but two or three times.
            %%% AND the second and third branch are not surface corrected!
            switch phase
               case {'PP','PKPPKP','PKIKPPKIKP'}
                   tt=tt+2*tt1;
               case {'PPP'}
                   tt=tt+4*tt1;
            end; % switch phase - inner switch for PP, PPP etc.
         else
            if strcmp(phase,'P')
               tt=tt2;
            else
               tt=NaN;
            end; % if strcmp
         end; % if takeoff<90
      end; % if isnan vdep
      %%% end of P, PP, ... %%%
      
   %%%% PKKP
   case {'PKKP'}
         if takeoff<=90
            %%% to arrive at the core, you have to go DOWN. takeoff<90!
            %%% compute vertex depth for K-leg
            vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
            if isnan(vdep)
               %%% NaN if no vertex exists
               %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Abort.']);
               tt=NaN;
               return;
            else
               %%% vertex exists.
               if isnan(model.cmb)
                    % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                    tt=NaN;
                    return;
               end; % if isnan(model.cmb)
               %%% find shallowest possible vertex below hypocenter, forget the other ones
               indies=find(vdep<=h); % was: (vdep<=model.rp-cmbspher); % below hypocenter & below CMB
               vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
               if vdep<r(length(r)-1)
                  %%% vertex is in deepest layer => abort.
                  %disp('MKX4P: vertex too deep. Abort.');
                  tt=NaN;
                  return;
               end;
               if isempty(vdep)
                  %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below ' num2str(cmbspher) 'km. Abort.']);
                  tt=NaN;
                  return;
               end; % if isempty(vdep)
               %%% FLAT EARTH transform for ray parameter
               p=mkpsfer2flat(p,rp);
               %%% FLAT EARTH transform for vdep and h (with dummy velocities)
               [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
               vdepflat=dmyz(1);
               hflat=dmyz(2);
               if (vdepflat<=model.cmb) %|(vdepflat>model.icb) % was: vdepflat<model.cmb
                  %%% vertex is outside core: this is not PKKP.
                  %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Abort.']);
                  tt=NaN;
                  return;
               end; % if vdep
               %%% compute P-leg traveltime for surface focus
               %disp('MKTIM4P: compute surface-surface time');
               tt1=mktpsum(p,vpflat,zflat,0,model.cmb,0);
               %%% compute P-leg traveltime from surface to real focus at HFLAT
               %disp('MKTIM4P: compute focal depth correction);
               tt2=mktpsum(p,vpflat,zflat,0,hflat,1);
               %%% compute K-leg trvaletime from CMB to vertex at VDEPFLAT
               %disp('MKTIM4P: compute K-leg traveltime');
               ttk=mktpsum(p,vpflat,zflat,model.cmb,vdepflat,0);
               %%% total traveltime
               %%% the following switch can be used to generate things like PKKKP etc as well.
               switch phase
                  case {'PKKP'}
                     %%% PKKP is one surface corrected P-leg, two K-legs and one CMB-to-surface
                     %%% P-leg.
                     tt=2*tt1-tt2+4*ttk;
               end; % inner switch
            end; % if isnan vdep
         else
            tt=NaN;
         end; % if takeoff<=90
      %%%% end of PKKP
      
      %%% PcP %%%%%%%%%%%%%%%%%
   case {'PcP','PcPPcP'}
      if takeoff<=90
         %%% you have to go DOWN to get reflected by the CMB!
          %% compute max allowed ray parameter
          %%% uses spherical model!!!
          %%% (bear in mind: discontinuity depths in MODEL are for flat, velocities are for spher!)
          maxp=mkraydepthinv(model.rp-cmbspher,model.vp,model.rp-model.z);
          if p>maxp
             %disp(['MKTIM4P: Ray parameter too large! p=' num2str(p) ', maxp=' num2str(maxp)]);
             tt=inf;
             return;
          end; % if p>maxp
          %%% compute vertex depth
          vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
          %%% if no vertex exists, stop computation and return inf.
          %%% But if one exists, it has to be below the CMB.
          %%% find shallowest possible vertex below hypocenter, forget the other ones
          indies=find(vdep<=h); % below or at hypocentral depth
          vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
          if isnan(model.cmb)
             %disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
             tt=NaN;
             return;
          end; % if isnan(model.cmb)
          if isempty(vdep)
             %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. returning inf.']);
             vdep=mkflat2sfer(model.cmb+1,model.cmb+1,model.rp); % dummy value, just to have something to transform
          end; % if isempty(vdep)
          %%% FLAT EARTH transform for ray parameter
          p=mkpsfer2flat(p,rp);
          %%% FLAT EARTH transform for vdep and h (with dummy velocities)
          [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
          vdepflat=dmyz(1);
          hflat=dmyz(2);
          if vdepflat<model.cmb
                %%% vertex is outside core: this is not PcP.
                %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not PcP. Returning inf.']);
                tt=inf;
                return;
          end; % if vdep
          %%% compute traveltime for surface focus
          tt1=mktpsum(p,vpflat,zflat,0,model.cmb,0);
          %%% compute traveltime from surface to real focus at HFLAT
          tt2=mktpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
          switch phase
             case {'PcP'}
                %%% total traveltime of phase is twice the time for surface focus
                %%% minus time from surface to focal depth
                tt=2*tt1-tt2;
             case {'PcPPcP'}
                %%% total PcPPcP traveltime is four times the time for surface focus
                %%% minus time from surface to focal depth
                tt=4*tt1-tt2;
          end; % inner switch
      else
         tt=NaN;
      end; % if takeoff<=90
   %%% end fo PcP %%%%%%%%%%
   
   
   %%% PcS, ScP %%%%%%%%%%%%%%%%%
   case {'PcS','ScP','PcSPcS','ScPScP'}
      if takeoff<=90
         %%% you have to go DOWN to get reflected by the CMB!
          %%% compute max allowed ray parameter
          %%% uses spherical model!!!
          %%% (bear in mind: discontinuity depths in MODEL are for flat, velocities are for spher!)
          maxp=mkraydepthinv(model.rp-cmbspher,model.vp,model.rp-model.z);
          if p>maxp
             %disp(['MKTIM4P: Ray parameter too large! p=' num2str(p) ', maxp=' num2str(maxp)]);
             tt=inf;
             return;
          end; % if p>maxp
          %%% compute vertex depth - it is sufficient ot have the PcP vertex below CMB
          vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
          %%% if no vertex exists, stop computation and return inf.
          %%% But if one exists, it has to be below the CMB.
          %%% find shallowest possible vertex below hypocenter, forget the other ones
          indies=find(vdep<=h); % below or at hypocentral depth
          vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
          if isnan(model.cmb)
             % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
             tt=NaN;
             return;
          end; % if isnan(model.cmb)
          if isempty(vdep)
             %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. returning inf.']);
             %tt=inf;
             %return; % set inf and stop computation. MK 03.07.2002
             vdep=mkflat2sfer(model.cmb+1,model.cmb+1,model.rp); % dummy value, just to have something to transform
          end; % if isempty(vdep)
          %%% FLAT EARTH transform for ray parameter
          p=mkpsfer2flat(p,rp);
          %%% FLAT EARTH transform for vdep and h (with dummy velocities)
          [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
          vdepflat=dmyz(1);
          hflat=dmyz(2);
          if vdepflat<model.cmb
                %%% vertex is outside core: this is not PcP.
                %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Returning inf.']);
                tt=inf;
                return;
          end; % if vdep
          %%% compute traveltime for surface focus
          ptt1=mktpsum(p,vpflat,zflat,0,model.cmb,0); % P-leg
          stt1=mktpsum(p,vsflat,zflat,0,model.cmb,0); % S-leg
          %%% compute traveltime from surface to real focus at HFLAT
          switch phase
             case {'PcS','PcSPcS'}
                tt2=mktpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
             case {'ScP','ScPScP'}
                tt2=mktpsum(p,vsflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
          end; % inner switch
          switch phase
             case {'PcS','ScP'}
                %%% total traveltime of phase is time for P-leg plus time for S-leg 
                %%% minus time from surface to focal depth
                tt=ptt1+stt1-tt2;
             case {'PcSPcS','ScPScP'}
                %%% total traveltime of phase is twice time for P-leg plus twice time for s-leg
                %%% minus time from surface to focal depth
                tt=2*ptt1+2*stt1-tt2;
          end; % inner switch
      else
         tt=NaN;
      end; % if takeoff<=90
   %%% end of PcS, ScP %%%%%%%%%%
   
   %%% PKiKP %%%%%%%%%%%%%%%%%
   case {'PKiKP'}
      if takeoff<=90
         %%% you have to go DOWN to rech the inner core!
          %%% compute vertex depth
          vdep=mkraydepth(p,vp,r,h,model.rp,disconradii); % vdep is a radius!
          %%% existence of vertices is more or less irrelevant for PcP. If no vertex exists - don't care.
          %%% But if one exists, it has to be below the CMB.
          %%% find shallowest possible vertex below hypocenter, forget the other ones
          indies=find(vdep<=h); % below or at hypocentral depth
          vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
          if isnan(model.cmb)|isnan(model.icb)
             % disp(['MKX4P: CMB or ICB not defined - no ' phase]);
             tt=NaN;
             return;
          end; % if isnan(model.cmb)
          if isempty(vdep)
             %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Setting dummy.']);
             vdep=mkflat2sfer(model.icb+1,model.icb+1,model.rp); % dummy value, just to have something to transform
          end; % if isempty(vdep)
          %%% FLAT EARTH transform for ray parameter
          p=mkpsfer2flat(p,rp);
          %%% FLAT EARTH transform for vdep and h (with dummy velocities)
          [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
          vdepflat=dmyz(1);
          hflat=dmyz(2);
          if model.icb-vdepflat>zprecision % MK13102006, was: vdepflat<model.icb
                %%% vertex is outside core: this is not PcP.
                %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not PcP. Abort.']);
                tt=inf;
                return;
          end; % if vdep
          %%% compute traveltime for surface focus
          tt1=mktpsum(p,vpflat,zflat,0,model.icb,0);
          %%% compute traveltime from surface to real focus at HFLAT
          tt2=mktpsum(p,vpflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
          %%% total traveltime of phase is twice the time for surface focus
          %%% minus time from surface to focal depth
          tt=2*tt1-tt2;
      else
         tt=NaN;
      end; % if takeoff<=90
   %%% end fo PKiKP %%%%%%%%%%
      
   case{'S','SS','SSS'} %%% S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% compute vertex depth
      vdep=mkraydepth(p,vs,r,h,model.rp,disconradii); % vdep is a radius!
      if isnan(vdep)&(takeoff<90) % takeoff<90: mk29092003
         %%% NaN if no vertex exists
         %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Abort.']);
         tt=NaN;
         return;
      else
         %%% vertex exists or ray is upgoing.
         %%% FLAT EARTH transform for ray parameter
         p=mkpsfer2flat(p,rp);
         %%% FLAT EARTH transform for h (with dummy velocities)
         [dmyv,dmyz]=mksfer2flat([1 1],rp-h,rp);
         hflat=dmyz;
         if takeoff<=90 % MK11102006, was: takeoff<90 % if takeoff<90 mk29092003
            %%% find shallowest possible vertex below hypocenter, forget the other ones
            indies=find(vdep<h); % below hypocenter
            vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
            if vdep<r(length(r)-1)
               %%% vertex is in deepest layer => abort.
               %disp('MKX4P: vertex too deep. Abort.');
               tt=inf;
               return;
            end;
            if isempty(vdep)
               %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg above ' num2str(h) 'km. Abort.']);
               tt=NaN;
               return;
            end; % if isempty(vdep)
            if vdep-(rp-cmbspher)<zprecision % was: vdep<rp-cmbspher
               %disp(['MKTIM4P: vertex inside core not allowed for ' phase '. Abort.']);
               tt=NaN;
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
                %%% S vertex is on the CMB and ray is not grazing!
                %%% This could become ScS if a P
                %%% vertext below the CMB exists.
                vdepp=mkraydepth(p,vp,r,h,model.rp,disconradii);
                indies=find(vdepp<=rp-cmbspher); 
                if ~isempty(indies)
                    %%% shallowest vertext below source is below CMB
                    %%% phase has a possible P vertex below the CMB -> not an S phase!
                    tt=NaN;
                    return;
                end; % if max(vdepp)>cmbspher
            end; % if vdep==rp-cmbspher
            %%%%% end of ScS suppression code MK28062006
            
         end; % if takeoff<90
         %%% compute traveltime from surface to real focus at HFLAT
         tt2=mktpsum(p,vsflat,zflat,0,hflat,1);
         if takeoff<=90 % MK11102006, was: takeoff<90
            %%% FLAT EARTH transform for vdep (with dummy velocities)
            [dmyv,dmyz]=mksfer2flat([1 1],rp-vdep,rp);
            vdepflat=dmyz;
            if vdepflat-model.cmb>zprecision % was: vdepflat>model.cmb MK30062006 % ">=" yields problems with S touching CMB. changed to ">" MK19052004
               %%% would penetrate the core, but the core is reserved for SKS and SKIKS
               %%% this is simply to avoid confusion, a calculation of S inside the core
               %%% would yield true paths and times!
               tt=NaN;
               return;
            end; % if vdepflat
            %%% compute traveltime for surface focus
            tt1=mktpsum(p,vsflat,zflat,0,vdepflat,0);
            %%% total traveltime of phase is twice the time for surface focus
            %%% minus time from surface to focal depth
            tt=2*tt1-tt2; 
            %%% generate SS and SSS if required
            %%% SS and SSS is the same as S - but two or three times.
            %%% AND the second and third branch are not surface corrected!
            switch phase
               case {'SS'}
                   tt=tt+2*tt1;
               case {'SSS'}
                   tt=tt+4*tt1;
            end; % switch phase - inner switch for PP, PPP etc.
         else
            if strcmp(phase,'S')
               tt=tt2;
            else
               tt=NaN;
            end; % if strcmp
         end; % if takeoff<90
      end; % if isnan vdep
      %%% end of S, SS, ... %%%%%
      
      %%% ScS %%%
   case {'ScS','ScSScS'}
      if takeoff<=90
         %%% you have to go DOWN to get reflected by the CMB!
         
          %%% catch some exceptional values
          if isnan(model.cmb)
             %disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
             tt=NaN;
             return;
          end; % if isnan(model.cmb)
         
          %%% compute max allowed ray parameter
          %%% uses spherical model!!!
          %%% (bear in mind: discontinuity depths in MODEL are for flat, velocities are for spher!)
          maxp=mkraydepthinv(model.rp-cmbspher,model.vs,model.rp-model.z);
          if p>maxp
             %%% ray misses the core (too shallow)
             tt=NaN;
             return;
          end; % if p>maxp
          
          
          %%% FLAT EARTH transform for ray parameter
          p=mkpsfer2flat(p,rp);
          %%% FLAT EARTH transform for vdep and h (with dummy velocities)
          [dmyv,hflat]=mksfer2flat([1 1],rp-[h],rp); 

          %%% compute traveltime for surface focus
          tt1=mktpsum(p,vsflat,zflat,0,model.cmb,0);
          %%% compute traveltime from surface to real focus at HFLAT
          tt2=mktpsum(p,vsflat,zflat,0,hflat,1); % novertex=1: suppress vertex code
          switch phase
             case {'ScS'}
                %%% total traveltime of ScS is twice the time for surface focus
                %%% minus time from surface to focal depth
                tt=2*tt1-tt2;
             case {'ScSScS'}
                %%% total traveltime of ScSScS is four times the time for surface focus
                %%% minus time from surface to focal depth
                tt=4*tt1-tt2;
          end; % inner switch
      else
         tt=NaN;
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
             %disp(['MKTIM4P: S leg vertex above CMB (' num2str(vdepsleg) ') - this is not a core phase. Returning inf.']);
             tt=inf;
             return;
          end; % if isempty(indies)
          %%% compute vertex depth for K-leg
          vdep=mkraydepth(p,vp,r,model.rp-cmbspher,model.rp,disconradii); % MK04102006 was: vdep=mkraydepth(p,vp,r,h,model.rp); % vdep is a radius!
          if isnan(vdep)
             %%% NaN if no vertex exists
             %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg. Abort.']);
             tt=NaN;
             return;
          else
             %%% vertex exists.
             %%% find shallowest possible vertex below CMN/ICB, forget the other ones
             switch phase
                case {'SKS','SKSSKS','SKKS'}
                   if isnan(model.cmb)
                     % disp(['MKX4P: Core Mantle Boundary not defined - no ' phase]);
                     tt=NaN;
                     return;
                   end; % if isnan(model.cmb)
                   indies=find(vdep<=model.rp-cmbspher); % below CMB
                case {'SKIKS'}
                   if isnan(model.cmb)|isnan(model.icb)
                     % disp(['MKX4P: CMB or ICB not defined - no ' phase]);
                     tt=NaN;
                     return;
                   end; % if isnan(model.cmb)
                   indies=find(vdep<=model.rp-icbspher); % below ICB
             end; % inner switch 
             vdep=max(vdep(indies)); % shallowest one (shallow == large radius)
             if vdep<r(length(r)-1)
                %%% vertex is in deepest layer => abort.
                %disp('MKX4P: vertex too deep. Abort.');
                tt=inf;
                return;
             end;
             if isempty(vdep)
                %disp(['MKTIM4P: no vertex for ' phase ' at p=' num2str(p*pi/180) 's/deg below ' num2str(cmbspher) 'km. Abort.']);
                tt=NaN;
                return;
             end; % if isempty(vdep)
             %%% FLAT EARTH transform for ray parameter
             p=mkpsfer2flat(p,rp);
             %%% FLAT EARTH transform for vdep and h (with dummy velocities)
             [dmyv,dmyz]=mksfer2flat([1 1],rp-[vdep h],rp);
             vdepflat=dmyz(1);
             hflat=dmyz(2);
             switch phase
                case {'SKS','SKSSKS','SKKS'}
                   if vdepflat<model.cmb % was: (vdepflat<=model.cmb) MK28092006, even before: %|(vdepflat>model.icb) % was:vdepflat<model.cmb
                      %%% vertex is outside core: this is not SKS.
                      %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of core (flat:' num2str(model.cmb) '): this is not ' phase '. Abort.']);
                      tt=inf;
                      return;
                   end; % if vdep
                case {'SKIKS'}
                   if (vdepflat<=model.icb) % was:vdepflat<model.cmb
                      %%% vertex is outside inner core: this is not SKIKS.
                      %disp(['MKTIM4P: vertex depth (flat) ' num2str(vdepflat) ' outside of inner core (flat:' num2str(model.icb) '): this is not ' phase '. Abort.']);
                      tt=inf;
                      return;
                   end; % if vdep
              end; % inner switch
             %%% compute S-leg traveltime for surface focus
             %disp('MKTIM4P: compute surface-surface time');
             tt1=mktpsum(p,vsflat,zflat,0,model.cmb,0);
             %%% compute S-leg traveltime from surface to real focus at HFLAT
             %disp('MKTIM4P: compute focal depth correction);
             tt2=mktpsum(p,vsflat,zflat,0,hflat,1);
             %%% compute K-leg trvaletime from CMB to vertex at VDEPFLAT
             %disp('MKTIM4P: compute K-leg traveltime');
             ttk=mktpsum(p,vpflat,zflat,model.cmb,vdepflat,0);
             %%% total traveltime
             switch phase
                case {'SKS','SKIKS'}
                   %%% total traveltime of phase is twice the time for surface focus
                   %%% minus time from surface to focal depth
                   %%% plus twice the k-leg traveltime from the core
                   %%% SKIKS is essentially the same as SKS
                   tt=2*tt1-tt2+2*ttk;
                case{'SKSSKS'}
                   %%% SKSSKS time is twice the time for surface-to-surface SKS
                   %%% minus surface correction
                   tt=2*(2*tt1+2*ttk)-tt2;
                case {'SKKS'}
                    %%% SKKS time is the same as SKS, but with an additional K-leg
                   tt=2*tt1-tt2+4*ttk;
             end; % inner switch: SKS, SKSSKS, SKKS
          end; % if isnan vdep
      else
         tt=NaN;
      end; % if takeoff<=90
      %%% end of SKS %%%
      
   otherwise %%% unkown phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      error(['MKTIM4P: cannot handle phase ' phase '. Abort.']);
      tt=NaN;
      return;
end; % switch phase
   

%toctime=toc; disp(['MKTIM4P: elapsed time: ' num2str(toctime) 's.']);

%%% return results
% tt=tt;
% vdep=vdep;