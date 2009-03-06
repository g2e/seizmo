function pscan=mkpsampler(phase,h,imodel);
% mkpsampler.....samples the Delta(p) function for given phase and source
%
% call: pscan=mkpsampler(phase,h,imodel);
%
%       phase: string containing a valid seismic phase name that can be
%       evaluated by MKX4P and MKTIM4p.
%
%           h: source depth below the surface [km]
%
%       imodel: MODEL structure as returned by MKREADND
%               OR
%               Improved MODEL structure as returned by MKIMPROVEMODEL
%
%              The criticalrays-field inserted by MKIMPROVEMODEL is of
%              utmost importance for the operation of MKPSAMPLER. However,
%              if it is not contained, it will be generated on the fly
%              (which means to waste CPU time!)
%
% result: pscan: [struct]
%                Strcutrue containing the complete scan result and the
%                positions at which continuous pieces of the delta(p)
%                function begin and end. The structure contains the
%                following fields:
%
%                .phase: [string]
%                        Name of the phase for which this struct is valid
%                .h: [km]
%                        focal depth for which this struct is valid
%                .angles: [deg]
%                         list of take off angles
%                .p: [sec/deg]
%                         list of ray parameters correpsonding to the take off
%                         angles at depth h
%                .dist: [deg]
%                         list of epicentral distances corresponding to p.
%                .vp: [km/s]
%                         P wave velocity at focal depth
%                .vs: [km/s]
%                         S wave velocity at focal depth
%                .starts: [index]
%                         positions where continuous pieces of dist(p) begin
%                .ends: [index]
%                         positions where continuous pieces of dist(p) end
%                         The i-th continuous piece begins at starts(i) and ends at
%                         ends(i).
%
% 
% This routine replaces the old MKSCANP routine.
% It produces a smarter sampling of the epicentral distance as function of
% ray parameter and take off angle as MKSCANP did. The main advantage is a
% more precise definition of the kinks and extrema of the distance
% function.
%
% Martin Knapmeyer, 23.06.2006

% 13122006 pscan entries with dist==inf are also removed

%%% init result
pscan=mkemptypscan(1);
pscan.phase=phase;
pscan.h=h;


%%% yet another magic number
%%% angleeps: used for comparison of resulting take off angles with angle
%%%           range boundaries found by MKSMARTTAKEOFF. For some reason
%%%           that I do not understand, the takeoff to the CMB differs by
%%%           some very small fraction of 1e-14, depending on the source.
%%%           So I blur the comparison by a somewhat larger epsilon.
angleeps=1e-6;


%%% does the IMODEL structure already contain the critical ray paramaters
%%% list? If not, create one!
if ~isfield(imodel,'criticalrays')
    %%% imodel does not contain a cirtical rays list, create one!
    imodel=mkimprovemodel(imodel);
end; % if ~isfield(imodel,'criticalrays')


%%% determine velocity at source depth (needed for angle-ray parm conversion)
%%% if h is a discontinuity depth, there will be two velocity values.
%%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
focus=mkinterpmodel(imodel,h,'simple');
indy=find(focus.z==h);
pscan.vp=focus.vp(indy);
pscan.vs=focus.vs(indy);


%%% Which critical ray parameters do we have to use?
%%% If phase name contains both P and S legs, both sets of ray parameters
%%% have to be used, if there is only one wave type, only a part needs to
%%% be used.
rayparmlist=[];
if (~isempty(strfind(lower(phase),'p'))) || (~isempty(strfind(lower(phase),'k')))
   %%% phase has P wave legs, include P ray parameters
   rayparmlist=[rayparmlist; imodel.criticalrays.p];
end; % if (~isempty(strfind(lower(phase),'p'))) | (~isempty(strfind(lower(phase),'k')))
if (~isempty(strfind(lower(phase),'s'))) || (~isempty(strfind(lower(phase),'j')))
   %%% phase has S wave legs, include S ray parameters
   rayparmlist=[rayparmlist; imodel.criticalrays.s];
end; %



%%% convert ray parameters to take off angles
takeofflist=mkrayp2angle(phase,h,rayparmlist,pscan.vp,pscan.vs,imodel.rp);
%dmyp=mkangle2rayp(phase,h,takeofflist,pscan.vp,pscan.vs,imodel.rp);



%%% the conversion of ray parameters into takeoff angles returns two angles
%%% per ray parameters, and the list of downgoing angles is concatenated at
%%% the end of the list of upgoing vectors.
%%% Therefore we double the ray parameters list to keep
%%% the two lists parallel and of identical size - but only if h~=0,
%%% because if h==0, only the downgoing rays are returned by MKRAYP2ANGLE.
if h~=0
   rayparmlist=[rayparmlist; rayparmlist];
end; % if h~=0



%%% remove all take off angles outside the angular interval possible for
%%% the current phase. also remove NaN angles.
[anglerange,minangle,maxangle]=mksmarttakeoff(phase,imodel,h,10);
%anglerange=0.1:90; minangle=0.1; maxangle=90; % test MK29092006
indies=find((takeofflist>=minangle-angleeps)&...
            (takeofflist<=maxangle+angleeps)&...
            (~isnan(takeofflist)));
%%% write results so far into pscan
pscan.angles=takeofflist(indies);
pscan.p=rayparmlist(indies);



%%% the min and max possible angle (acc. to MKSMARTTAKEOFF) is typically not found
%%% in the critical rays analysis. Therefore we append them by hand here.
pscan.angles=[pscan.angles; minangle; maxangle];
pscan.p=[pscan.p; mkangle2rayp(pscan.phase,pscan.h,[minangle; maxangle],imodel)];
[newminp,minangle]=mkfumblep(pscan.p(end-1),pscan.angles(end-1),'both',pscan.phase,pscan.h,imodel);
[newmaxp,maxangle]=mkfumblep(pscan.p(end),pscan.angles(end),'both',pscan.phase,pscan.h,imodel);
pscan.p=[pscan.p; newminp; newmaxp];
pscan.angles=[pscan.angles; minangle; maxangle];
% was: (MK10102006)
% pscan.p(end-1)=newminp;
% pscan.angles(end-1)=minangle;
% pscan.p(end)=newmaxp;
% pscan.angles(end)=maxangle;


%%% another important takeoff angle is 90 degrees, especially since poor
%%% depth sampling of the velocoty model causes strange things around
%%% there.
pscan.angles=[pscan.angles; 90];
pscan.p=[pscan.p; mkangle2rayp(pscan.phase,pscan.h,90,imodel)];


%%% scan over all critical takeoff angles from model analysis 
takeoffanz=length(pscan.angles);
pscan.dist=pscan.angles*0;
%pscan.maxz=pscan.dist;
for indy=1:takeoffanz
    dist=mkx4p(pscan.phase,pscan.h,pscan.p(indy),imodel,pscan.angles(indy));
    pscan.dist(indy)=dist;
    %pscan.maxz(indy)=max(segz);
    %tt(indy)=mktim4p(pscan.phase,pscan.h,pscan.p(indy),imodel,pscan.angles(indy));
end; % for indy


%%% sort entries by takeoff angle
%%% this is necessacry because the minimum/maximum search below assumes
%%% that the list is ordered.
[pscan.angles,sorter]=sort(pscan.angles);
pscan.p=pscan.p(sorter);
pscan.dist=pscan.dist(sorter);


% %%% control plot
% subplot(1,2,1);
% plot(pscan.angles,pscan.dist,'bo-');


%%% test for the existence of NaN distance values in PSCAN. Such values
%%% make the later interpolation, er, difficult.
pscan=mknanfreepscan(pscan);


%%% search for additional local minima or maxima where the first derivative
%%% of delta(alpha) becomes zero
%%% Each pair of adjacent samples in the current PSCAN struct defines an
%%% interval that might contain a local minimum or maximum. These intervals
%%% are thus tested subsequently.
intervals=length(pscan.angles)-1;
for intervalcnt=1:intervals
    
    %%% define current interval (column vector!!)
    alphalimit=[pscan.angles(intervalcnt); pscan.angles(intervalcnt+1)];
    deltalimit=[pscan.dist(intervalcnt); pscan.dist(intervalcnt+1)];
    
%     %% control plot
%     subplot(1,2,1);
%     %plot(pscan.angles,pscan.dist,'bo-');
%     hold on
%     plot(alphalimit,deltalimit,'rs');
%     hold off
    
    %%% search for local extremum in the current interval
    [newalpha,newp,newdelta]=mkfindextremedists(alphalimit,deltalimit,phase,h,imodel);
    
    %%% append solution to PSCAN struct
    pscan.angles=[pscan.angles; newalpha];
    pscan.p=[pscan.p; newp];
    pscan.dist=[pscan.dist; newdelta];
    
end; % for intervalcnt


%%% The above loop may have generated new NaN distances. therefore...
%%% test for the existence of NaN distance values in PSCAN. Such values
%%% make the later interpolation, er, difficult.
pscan=mknanfreepscan(pscan);



%%% sort entries by takeoff angle
%%% this is necessacry because the results of the search above were just
%%% appended unordered.
[pscan.angles,sorter]=sort(pscan.angles);
pscan.p=pscan.p(sorter);
pscan.dist=pscan.dist(sorter);

%%% remove repeated samples
[pscan.angles,pscan.p,pscan.dist]=mkuniquesamples(pscan.angles,pscan.p,pscan.dist);

%%% construct starts and ends lists
%%% these are essentially historical artefacts from the MKSCANP shooting
%%% and constructed only to provide backward compatibility of the output.
pscan.starts=1:(length(pscan.angles)-1);
pscan.ends=2:length(pscan.angles);

return;


