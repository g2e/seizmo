function newmod=mkinterpmodel(model,newz,preserve);
% mkinterpmodel..........depth interpolation of velocity model
%
% call: newmod=mkinterpmodel(model);
%       newmod=mkinterpmodel(model,newz);
%       newmod=mkinterpmodel(model,newz,preserve);
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
%         newz:   new z values to which MODEL is interpolated
%                 This is assumed to be a vector ot type min:delta:max with no repetitions,
%                 or a scalar value!
%                 CAUTION: it is important to leave all 1st order discontinuities of the model intact.
%                          MKINTERPMODEL therefore adds all discontinuity depths to the
%                          z value list in NEWZ to guarantee that discontinuities remain there.
%                 If NEWZ is not present, interpolated layers will have a constant thickness of 10km
%                 in flat coordinates, except for the lowermost layer, which contains the
%                 planetary centera d goes down to inf. Enter the empty matrix [] if you wish to get this
%                 automatic interpolation in three-argument-mode.
%         preserve: using this switch, you can disable the mentioned discontinuity preservation.
%                   possible values:
%                   'preserve': to leave discontinuities intact
%                   'simple': for a simple interpolation which does not care for discontinuities
%                             use this when preparing a model file for PSM modeling or when you
%                             want to know a single velocity or densitiy value at a single depth.
%                             Do NOT use this if you wish to compute rays of travel times!
%                   DEFAULT: 'preserve' (used in two-arguments version)
%
% result: newmod: a new model structure in which newmod.z=NEWZ, and vp,vs,rho,qp and qs are
%                 interpolated to NEWZ.
%
% It is difficult to identify higher-order discontinuities by looking at the model.z list.
% MKINTERPMODEL explicitly preserves the follwing types of discontinuity:
% - all discontinuities defined by the .conr, .moho etc fields
% - all discontinuities listed in dz and dname
% - all unnamed discontinuities that are not mentioned in dz and dname, but at which
%   double velocity definitions are given (two veolocities at the same depths, i.e. two adjacent
%   depth values in model.z are identical)
%
% A fine sampling of the model is important mainly when you deal with small ray parameters
% and rays that nearly reach the planet's center: this is somewhere near infinite depth
% in Flat Earth Transformation and thus larger numerical errors occur. I made good experiences
% with a sampling that produces 5km layers under the FET. But for IASP91, for example, such
% dense sampling is necessary only in the inner core. You don't need to produce a dense
% sampling at all depths - and that's good, because densely sampled models require more
% computation time in ray path and travel time evaluation.
%
% Martin Knapmeyer, 28.05.2002, 05.07.2002, 18.07.2002
%                   completely rewritten 08.09.2003
%                   16.06.2006
%
% BUGS: may produce nonsense when you try to interpolate to a list
%       containing repeating depths (e.g. 10 20 20 30). (This does not
%       affect existing repetitions from discontinuities, but you cannot
%       introduce new repetitions this way)



%%% understand input
if (nargin==1)|(isempty(newz))
   preserve='preserve';
   %%% automatic interpolation: construct a z-vector which is equal spaced in
	%%% flat coordinates
	dz=5; % [km] - spacing of equidistant layers in flat coordinates
	disp(['MKINTERPMODEL: autospacing: ' num2str(dz) 'km layer thickness in flat coordinates']);
	[dmy,rpf]=mksfer2flat(model.rp,model.z(length(model.z))-1,model.rp); %mksfer2flat(model.rp,model.z(length(model.z)-1),model.rp);
	zf=0:dz:rpf;
	[dmy,zf]=mkflat2sfer(zf,zf,model.rp); % back into sphercial coordinates
	newz=zf;
end; % if nargin
if nargin==2
   preserve='preserve';
end; % if nargin
preserve=lower(preserve);


%%% NEWZ is the vector of depths at which the interpolated model has to be
%%% defined.
%%% make NEWZ a col vector - just to be sure
%%% and sort it! MK16062006
newz=sort(newz(:));

%%% evaluate PRESERVE
switch preserve
    case {'simple'}
         %%% nothing special happens
         %%% except for a warning MK07122006 (removed 11122006)
         %%% disp('MKINTERPMODEL: Models interpolated in SIMPLE mode must not be used for travel time or ray computations!');
    case {'preserve'}
    
         %%% find all first order discontinuities 
         %%% even the unnamed ones which are not listed in dz & dname
         %%% DISCINDY will contain index of elements pointing to upper-side data
         discindy=find(diff(model.z)==0);
         discz=model.z(discindy); % depth of all dubly defined depths
         %%% depths of all named and unnamed discontinuities
         discz=[discz(:)',...
                model.conr,...
                model.moho,...
                model.d410,...
                model.d520,...
                model.d660,...
                model.cmb,...
                model.icb,...
                model.dz(:)'];
         %%% remove NaNs
         discz=discz(find(~isnan(discz)));
         %%% make list unique
         discz=unique(sort(discz));
         
         %%% insert discontinuity depths list into NEWZ
         %%% but without generatin threefold or quadruple (or n-tuple) entries!
         %%% to do so, remove all entries from NEWZ which correspond to
         %%% entries in DISCZ
         for indy=1:length(discz)
             indies=find(newz==discz(indy)); % identify entries to be removed
             newz(indies)=newz(indies)+NaN;  % and mark them
         end; % for indy
         newz=newz(find(~isnan(newz)));  % this is the entry removal
         
         %%% generate double entry for each depth
         %%% this is needed since we need two lines with identical z
         %%% in the final model.z array
         %disczdoubled=sort([discz discz]); % not necessary MK26062006
         
         %%% insert into NEWZ
         %newz=sort([newz(:)' disczdoubled(:)'])'; % not necessary MK26062006
         newz=sort([newz(:)' discz(:)'])';
         
    otherwise
        error(['MKINTERPMODEL: unknown option ' upper(preserve)]);
end; % switch preserve


%%%%%% start interpolation
%%% NEWZ is now the vector of z values for which we wish to interpolate
%%% we now identify all pieces of model.z in which there are no repetitions
%%% of z values.

%%% but first, some simple things
newmod.z=[]; %newz;
newmod.vp=[]; %newz*0;
newmod.vs=[]; %newz*0;
newmod.rho=[]; %newz*0;
newmod.qp=[]; %newz*0;
newmod.qs=[]; %newz*0;
newmod.conr=model.conr;
newmod.moho=model.moho;
newmod.d410=model.d410;
newmod.d520=model.d520;
newmod.d660=model.d660;
newmod.cmb=model.cmb;
newmod.icb=model.icb;
newmod.dz=model.dz;
newmod.dname=model.dname;
newmod.year=model.year;
newmod.rp=model.rp;
newmod.name=model.name;


%%% the interpolation
%%% partition MODEL into pieces with unique z values
%%% these pieces are delimited by the discontinuities
indies=find(diff(model.z)==0); % top sides of all discontinuities in MODEL
modstartindy=[1; indies+1]';
modendindy=[indies; length(model.z)]';
pieceanz=length(modstartindy); % number of pieces in MODEL
%%% find corresponding pieces in NEWMOD
indies=find(diff(newmod.z)==0);
newstartindy=[1; indies+1]';
newendindy=[indies; length(newmod.z)]';
newpieceanz=length(newstartindy); % number of pieces in NEWMOD

%%% loop over all pieces that can interpolated at once.
%%% Do not use "spline" interpolation, since it would produce artefact LVZs
%%% around corners of the original model. "pchip" avoids that. 
method='pchip';
for piececnt=1:pieceanz
   
    %%% range of MODEL used for interpolation
    modrange=modstartindy(piececnt):modendindy(piececnt);
    %%% range of NEWZ for which we interpolate
    if newpieceanz>1
%        newrange=newstartindy(piececnt):newendindy(piececnt); %% old code:  removed 16062006
       
       %%% here follows new code MK16062006
       
       %%% now identify the interval of newmod.z that belongs to the
       %%% current depth interval defined by MODERANGE
       %%% we assume that if NEWZ contains repeating depths, these
       %%% constitute the discontinuities already present in MODEL, i.e.
       %%% it is not possible to inroduce new discontinuities via NEWZ.
       %%% But we also assume that NEWZ does not contain all
       %%% discontinuities present in MODEL.
       
       %%% find first NEWZ interval corresponding to current MODEL interval
       indies1=find(...
                   newmod.z(newstartindy)==model.z(modstartindy(piececnt))...
                  );
       if ~isempty(indies1)
           
           %%% there are intervals that start with the current depth
           %%% intervals starting depth
           
           firstnewsample=newstartindy(indies1(1));

           %%% find last sample in NEWZ that is inside current depth interval
           %%% AND inside the NEWZ interval beginning at NEWSTARTINDY
           subspace=(firstnewsample+1):length(newmod.z); %(firstnewsample+1):newendindy(indies1(1));
           indies2=find(...
                       newmod.z(subspace)<=model.z(modendindy(piececnt))...
                      );
           
           
           if ~isempty(indies2)
               %%% There are more sample belonging to the current depth
               %%% interval
               lastnewsample=subspace(indies2(end));
               
           else
               
               %%% no more samples belonging to current depth interval
               %%% first sample is also the last sample.
               lastnewsample=firstnewsample;
               
           end; % if ~isempty(indies2)

           %%% index range pointing to NEWZ, designating samples belonging to
           %%% the current depth interval
           newrange=firstnewsample:lastnewsample;
           
       else
           
           %%% There is no interval that starts with the starting depth of
           %%% the current depth interval
           %%% set NEWRANGE empty to suppress the interpolation below.
           newrange=[];
       
       end; % if ~isempty(indies1)
       
       %%% end of new code MK16062006
       
       
    else
       %%% NEWZ does not consist of distinct pieces - we have to fish out
       %%% the samples belonging to the current interval individually.
       newrange=find((newz>=min(model.z(modrange)))&(newz<=max(model.z(modrange))));
    end; % if newpieceanz
    
    
%     %%% control plot: what intervals have been selected?
%     figure(2);
%     clf;
%     plot(model.vp,model.z,'k');
%     hold on
%     plot(model.vs,model.z,'k');
%     plot(newmod.vp,newmod.z,'r.');
%     plot(newmod.vs,newmod.z,'r.');
%     plot(model.vp(modrange),model.z(modrange),'g.-');
%     plot(newmod.vp(newrange)+10,newmod.z(newrange),'bo-');
%     plot(model.vs(modrange),model.z(modrange),'g.-');
%     plot(newmod.vs(newrange)+5,newmod.z(newrange),'bo-');
%     hold off
%     drawnow;
    
      
    
    %%% the interpolation itself
    %%% interpolation is carried out only if newrange is not empty, because
    %%% only then the corresponding model piece overlaps with the desired z
    if ~isempty(newrange) % was: ~isnan(newrange) changed MK16062006
       %%% only those fields are sent to interp1 which are not completly NaN.
       %%% completely NaN fields (as rho, qp, qs for IASP91) are just extended
       %%% to longer completely NaN fields. MatLabR14 update, MK24.08.2004
       %%% condition to test: prod(real(isnan(SomeVector))) is zero if there are
       %%% non-NaN elements
       
       %%% z
       newmod.z=[newmod.z; newz(newrange)];

       %%% vp
       if prod(real(isnan(model.vp(modrange))))==0
          %%% not everything is NaN: interpolate to new depths
          interpresult=interp1(model.z(modrange),...
                                      model.vp(modrange),...
                                      newz(newrange),...
                                      method);
          newmod.vp=[newmod.vp; interpresult];
       else
          %%% everything is NaN: replace by chain of NaNs of new length
          newmod.vp=[newmod.vp; newz(newrange)*NaN];
       end; % if sum(isnan())

       %%% vs
       if prod(real(isnan(model.vs(modrange))))==0
          %%% not everything is NaN: interpolate to new depths
          interpresult=interp1(model.z(modrange),...
                                   model.vs(modrange),...
                                   newz(newrange),...
                                   method);
          newmod.vs=[newmod.vs; interpresult]; 
       else
          %%% everything is NaN: replace by chain of NaNs of new length
          newmod.vs=[newmod.vs; newz(newrange)*NaN];
       end; % if sum(isnan())   

       %%% rho
       if prod(real(isnan(model.rho(modrange))))==0
          %%% not everything is NaN: interpolate to new depths
          interpresult=interp1(model.z(modrange),...
                                   model.rho(modrange),...
                                   newz(newrange),...
                                   method);
          newmod.rho=[newmod.rho; interpresult];
       else
          %%% everything is NaN: replace by chain of NaNs of new length
          newmod.rho=[newmod.rho; newz(newrange)*NaN];
       end; % if sum(isnan())   

       %%% qp
       if prod(real(isnan(model.qp(modrange))))==0
          %%% not everything is NaN: interpolate to new depths
          interpresult=interp1(model.z(modrange),...
                                   model.qp(modrange),...
                                   newz(newrange),...
                                   method);
          newmod.qp=[newmod.qp; interpresult]    ;             
       else
          %%% everything is NaN: replace by chain of NaNs of new length
          newmod.qp=[newmod.qp; newz(newrange)*NaN];
       end; % if sum(isnan())  

       %%% qs
       if prod(real(isnan(model.qs(modrange))))==0
          %%% not everything is NaN: interpolate to new depths
          interpresult=interp1(model.z(modrange),...
                                   model.qs(modrange),...
                                   newz(newrange),...
                                   method);
          newmod.qs=[newmod.qs; interpresult];
       else
          %%% everything is NaN: replace by chain of NaNs of new length
          newmod.qs=[newmod.qs; newz(newrange)*NaN];
       end; % if sum(isnan())
    
    end; % if ~isnan(newrange)
    
    
   
   
end; % for piececnt


%%% make result col vectors
newmod.z=newmod.z(:);
newmod.vp=newmod.vp(:);
newmod.vs=newmod.vs(:);
newmod.rho=newmod.rho(:);
newmod.qp=newmod.qp(:);
newmod.qs=newmod.qs(:);


return;