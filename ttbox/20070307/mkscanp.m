function pscan=mkscanp(phase,h,model,dangle);
% MKSCANP........computes distance as function of ray parm
%
% call: pscan=mkscanp(phase,h,model);
%       pscan=mkscanp(phase,h,model,dangle);
%
%       phase: string containing seismic phase name like 'P', 'S', 'ScS', 'PKPdf', etc.\
%              Phase names are case sensitive!
%           h: focal depth [km]
%          model: A structure describing the velocity distribution.
%                 Such a structure can be obtained via MKREADND.
%      dangle: delta angle parameter: This gives the take off angle
%              increment for the scan, in degrees. Should be small.
%              default: 0.1;
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
%                    focal depth for which this struct is valid
%                .angles: [deg]
%                         list of take off angles
%                .p: [sec/deg]
%                    list of ray parameters correpsonding to the take off
%                    angles at depth h
%                .dist: [deg]
%                    list of epicentral distances corresponding to p.
%                .vp: [km/s]
%                    P wave velocity at focal depth
%                .vs: [km/s]
%                    S wave velocity at focal depth
%                .starts: [index]
%                    positions where continuous pieces of dist(p) begin
%                .ends: [index]
%                    positions where continuous pieces of dist(p) end
%                    The i-th continuous piece begins at starts(i) and ends at
%                    ends(i).
%               
%
% A ray which starts at depth H with a take of angle of ANGLES(i) will have
% ray parameter p(i) and arrive at an epicentral distance dist(i).
%
% This was a part of MKFINDP. But sometimes it is quicker to do this scan
% outside of MKFINDP, for later re-use.
%
% Martin Knapmeyer, 27.05.2005, 30.50.2005, 26.10.2006


%%% prepare table of delta as function of p
if nargin<4
    dangle=0.1; % angle increment. MUST BE THAT SMALL - OTHERWISE YOU'LL MISS SOLUTIONS!
                % Make this smaller (0.001) to get the now missing solutions
                % for P at 2dg distance in IASP91 or for the missing PKIKP at
                % 114-120deg. MK26102004
end; % if nargin
angles=mksmarttakeoff(phase,model,h,dangle); % angles at which source radiates


%%% determine velocity at source depth (needed for angle-ray parm conversion)
%%% if h is a discontinuity depth, there will be two velocioty values.
%%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
focus=mkinterpmodel(model,h,'simple');
indy=find(focus.z==h);
vp=focus.vp(indy);
vs=focus.vs(indy);

%%% convert angles into ray parameters
anz=length(angles);
dist=zeros(1,anz);
p=dist;
maxz=zeros(size(dist));
for indy=1:anz
   %%% at a discontinuity, there will be two velocity values.
   %%% use above-value for angle>=90 and below-value for angle<90   
   if angles(indy)>=90
      whichangle=1;
   else
      whichangle=min(2,length(vp));
   end; % if angle
   p(indy)=mkangle2rayp(phase,h,angles(indy),vp(whichangle),vs(whichangle),model.rp);
   %%% now we have ray parameters to compute a distance from
   [dist(indy),rx,rz]=mkx4p(phase,h,p(indy),model,angles(indy));
   maxz(indy)=max(rz);
end; % for angle
%%% now we have calculated p as function of angle and dist as function of p.


% %%% scan result as function of penetration depth
% subplot(1,2,2);
% cla;
% plot(dist,maxz,'.-');
% ylabel('Penetration Depth [km]');
% xlabel('Epicentral Distance [deg]');
% drawnow;


%%% the next task is to find any value of p for which dist(p) has the value
%%% given in input parameter DELTA. That's hard, since delta(p) is not smooth,
%%% not monotonic and not continuous.
%%% we have to cut delta(p) into continuous pieces and handle each of these pieces
%%% separately.
allreals=find((~isinf(dist))&(~isnan(dist)));
if ~isempty(allreals)
   lastreal=allreals(end);
	firstreal=allreals(1); % index of the first non-NaN-non-inf-Element of dist
	if firstreal==1
	   %starts=[firstreal find(diff((~isnan(dist))&(~isinf(dist)))>0)+1]; % indexes where continuous pieces begin
       starts=[firstreal find(diff((~isnan(dist))+(~isinf(dist)))>0)+1]; % new MK08092004
	else
	   %starts=[find(diff((~isnan(dist))&(~isinf(dist)))>0)+1];
       starts=[find(diff((~isnan(dist))+(~isinf(dist)))>0)+1]; % new MK08092004
	end; % if firstreal==1
	%ends=find(diff((isnan(dist))&(isinf(dist)))>0); % indexes where continuous pieces end
    ends=find(diff((isnan(dist))+(isinf(dist)))>0); % new MK08092004
	if length(ends)<length(starts)
       %%% the above procedure to find ends does not find the last element of DIST as end.
       %%% if the last continuous piece does not end before the last element, we have to 
       %%% add that to the list.
       ends=[ends lastreal]; %length(dist)];
	end; % if length(ends)
    
    %%% plot scan result as function of takeoff angle
	%clf; plot(angles,dist,'.-'); hold on; plot(angles(starts),dist(starts),'r>'); plot(angles(ends),dist(ends),'r<'); hold off; grid on
	%%% plot scan result as function of ray parameter
	%clf; plot(p,dist,'.-'); hold on; plot(p(starts),dist(starts),'r>'); plot(p(ends),dist(ends),'r<'); hold off; grid on
	
	
	
	%%% MK28102003
	%%% no we have identified all holes of function dist(p). But there might be
	%%% jumps in the function, which are not identified by the above procedure.
	%%% We therefore have to search for jumps in each of the pieces identified
	%%% up to now. Their begin and end indexes will be added to STARTS and ENDS.
	anz=length(starts);
	newstarts=[];
	newends=[];
	for i=1:anz
	    %%% mkfindjump finds the index of the last element before a jump
	    indexrange=starts(i):ends(i);
	    if ~isempty(indexrange)
	       jumps=mkfindjump(angles(indexrange),dist(indexrange));
	       jumps=jumps+(indexrange(1)-1); % make it an index into the complete array
	       additionalstarts=jumps+1;
	       additionalends=jumps;
	       newstarts=[newstarts starts(i) additionalstarts additionalends];
	       newends=[newends ends(i) additionalends additionalstarts];
	    end; % if ~isempty(indexrange)
	end; % for i
	starts=newstarts;
	ends=newends;
	starts=sort(starts);
	ends=sort(ends);
	
	%%% plot scan result as function of takeoff angle
	%clf; plot(angles,dist,'.-'); hold on; plot(angles(starts),dist(starts),'r>'); plot(angles(ends),dist(ends),'r<'); hold off; grid on
	%%% plot scan result as function of ray parameter
	%clf; plot(p,dist,'.-'); hold on; plot(p(starts),dist(starts),'r>'); plot(p(ends),dist(ends),'r<'); hold off; grid on

else
   %%% no arrivals exist!
   %p=NaN; %[]; % commented out MK26102006
   a=NaN;
   d=NaN; %[];
   starts=NaN; %MK18052005
   ends=NaN; %MK18052005
   deltain=NaN; %MK18052005
end; % if ~isempty(allreals)

%%% built result struct
pscan=mkemptypscan(1);
pscan.phase=phase;
pscan.h=h;
pscan.angles=angles;
pscan.p=p;
pscan.dist=dist;
pscan.vp=vp;
pscan.vs=vs;
pscan.starts=starts;
pscan.ends=ends;

