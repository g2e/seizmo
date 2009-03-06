function [p,a,d,deltain,starts,ends]=mkfindp(phase,delta,h,model,pscan);
% mkfindp.........finds ray parameter corrsponding to an epicentral distance IN SPHERICAL EARTH
%
% call: [p,a,d,deltain,starts,ends]=mkfindp(phase,delta,h,model);
%       [p,a,d,deltain,starts,ends]=mkfindp(phase,delta,h,model,pscan);
%
%       phase: string containing seismic phase name like 'P', 'S', 'ScS', 'PKPdf', etc.\
%              Phase names are case sensitive!
%       delta: epicentral distance [deg]
%              might be a vector of epicentral distances
%           h: focal depth [km]
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
%          pscan: output of mkscanp(pahse,h,model);
%                 If these parameters are not given, MKSCANP will be called
%                 internally. External computation sometimes saves time.
%
%
% result: p: ray parameter [sec/deg] at which the given PHASE appears at distance DELTA
%            NaN if no such ray parameter exists (shadow zones etc)
%         a: take off angle of the PHASE at the source
%            necessary to distinguish between upgoing and downgoing rays
%         d: distances at which the respective rays arrive [deg]
%            you need this since the ray parameter list might be longer than DELTA list.
%            (and there might be small differences between the input values and these
%             results)
%         deltain: list of delta values for which the result is valid
%                  p(i), a(i), d(i) are the results for input value deltain(i)
%         starts: all index values (index to p,a,d) at which a continuous piece
%                 of function d(a) begins
%         ends:   all index values (index to p,a,d) at which a continuous piece
%                 of function d(a) ends.
%
% iterative search procedure.
%
% Martin Knapmeyer, 24.04.2002, 05.07.2002, 23.09.2003, 25.09.2003
%                   return NaN if no solution: 14.10.2003
%                   return start, ends, improved cutting: 28.10.2003
%                   scan phase moved to MKSCANP 27.05.2005, 30.05.2005

%tic;

error('MKFINDP: the use of MKFINDP is deprecated, please use MKSHOOTRAY instead!');

%%% initialize result
pres=[];
ares=[];
dres=[];
a=[];
d=[];
deltain=[];


%%% some constants
radian=pi/180;
epsilondeg=0.01; % epsilon used in delta-search [deg]

%%% call MKSCANP if necessary
if nargin==4
    p=NaN;
    pscan=mkscanp(phase,h,model);
end; % if nargin==4

%%% decompose pscan (necessary because of history of this function...
angles=pscan.angles;
p=pscan.p;
dist=pscan.dist;
vp=pscan.vp;
vs=pscan.vs;
starts=pscan.starts;
ends=pscan.ends;


%%% options for MKFINDZEROS: control parameters
fzopt.initwidth=0.015;
fzopt.wdtdivisor=1/(0.5*(sqrt(5)-1));
fzopt.epsilon=0.001;
fzopt.maxcnt=50;
%disp(['MKFINDP: initwidth: ' num2str(fzopt.initwidth) ', wdtdivisor: ' num2str(fzopt.wdtdivisor)]);

% %%% prepare table of delta as function of p
% dangle=0.1; % angle increment. MUST BE THAT SMALL - OTHERWISE YOU'LL MISS SOLUTIONS!
%             % Make this smaller (0.001) to get the now missing solutions
%             % for P at 2dg distance in IASP91 or for the missing PKIKP at
%             % 114-120deg. MK26102004
% angles=mksmarttakeoff(phase,model,h,dangle); % angles at which source radiates
% 
% 
% %%% determine velocity at source depth (needed for angle-ray parm conversion)
% %%% if h is a discontinuity depth, there will be two velocioty values.
% %%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
% focus=mkinterpmodel(model,h,'simple');
% indy=find(focus.z==h);
% vp=focus.vp(indy);
% vs=focus.vs(indy);
% 
% %%% convert angles into ray parameters
% anz=length(angles);
% dist=zeros(1,anz);
% p=dist;
% for indy=1:anz
%    %%% at a discontinuity, there will be two velocity values.
%    %%% use above-value for angle>=90 and below-value for angle<90   
%    if angles(indy)>=90
%       whichangle=1;
%    else
%       whichangle=min(2,length(vp));
%    end; % if angle
%    p(indy)=mkangle2rayp(phase,h,angles(indy),vp(whichangle),vs(whichangle),model.rp);
%    %%% now we have ray parameters to compute a distance from
%    [dist(indy),rx,rz]=mkx4p(phase,h,p(indy),model,angles(indy));
% end; % for angle
% %%% now we have calculated p as function of angle and dist as function of p.

%toc

% %%% the next task is to find any value of p for which dist(p) has the value
% %%% given in input parameter DELTA. That's hard, since delta(p) is not smooth,
% %%% not monotonic and not continuous.
% %%% we have to cut delta(p) into continuous pieces and handle each of these pieces
% %%% separately.
allreals=find((~isinf(dist))&(~isnan(dist)));
if ~isempty(allreals)
%    lastreal=allreals(end);
% 	firstreal=allreals(1); % index of the first non-NaN-non-inf-Element of dist
% 	if firstreal==1
% 	   %starts=[firstreal find(diff((~isnan(dist))&(~isinf(dist)))>0)+1]; % indexes where continuous pieces begin
%        starts=[firstreal find(diff((~isnan(dist))+(~isinf(dist)))>0)+1]; % new MK08092004
% 	else
% 	   %starts=[find(diff((~isnan(dist))&(~isinf(dist)))>0)+1];
%        starts=[find(diff((~isnan(dist))+(~isinf(dist)))>0)+1]; % new MK08092004
% 	end; % if firstreal==1
% 	%ends=find(diff((isnan(dist))&(isinf(dist)))>0); % indexes where continuous pieces end
%     ends=find(diff((isnan(dist))+(isinf(dist)))>0); % new MK08092004
% 	if length(ends)<length(starts)
%        %%% the above procedure to find ends does not find the last element of DIST as end.
%        %%% if the last continuous piece does not end before the last element, we have to 
%        %%% add that to the list.
%        ends=[ends lastreal]; %length(dist)];
% 	end; % if length(ends)
%     
%     %%% plot scan result as function of takeoff angle
% 	%clf; plot(angles,dist,'.-'); hold on; plot(angles(starts),dist(starts),'r>'); plot(angles(ends),dist(ends),'r<'); hold off; grid on
% 	%%% plot scan result as function of ray parameter
% 	%clf; plot(p,dist,'.-'); hold on; plot(p(starts),dist(starts),'r>'); plot(p(ends),dist(ends),'r<'); hold off; grid on
% 	
% 	
% 	
% 	%%% MK28102003
% 	%%% no we have identified all holes of function dist(p). But there might be
% 	%%% jumps in the function, which are not identified by the above procedure.
% 	%%% We therefore have to search for jumps in each of the pieces identified
% 	%%% up to now. Their begin and end indexes will be added to STARTS and ENDS.
% 	anz=length(starts);
% 	newstarts=[];
% 	newends=[];
% 	for i=1:anz
% 	    %%% mkfindjump finds the index of the last element before a jump
% 	    indexrange=starts(i):ends(i);
% 	    if ~isempty(indexrange)
% 	       jumps=mkfindjump(angles(indexrange),dist(indexrange));
% 	       jumps=jumps+(indexrange(1)-1); % make it an index into the complete array
% 	       additionalstarts=jumps+1;
% 	       additionalends=jumps;
% 	       newstarts=[newstarts starts(i) additionalstarts additionalends];
% 	       newends=[newends ends(i) additionalends additionalstarts];
% 	    end; % if ~isempty(indexrange)
% 	end; % for i
% 	starts=newstarts;
% 	ends=newends;
% 	starts=sort(starts);
% 	ends=sort(ends);
% 	
% 	%%% plot scan result as function of takeoff angle
%	clf; plot(angles,dist,'.-'); hold on; plot(angles(starts),dist(starts),'r>'); plot(angles(ends),dist(ends),'r<'); hold off; grid on
% 	%%% plot scan result as function of ray parameter
% 	clf; plot(p,dist,'.-'); hold on; plot(p(starts),dist(starts),'r>'); plot(p(ends),dist(ends),'r<'); hold off; grid on
	
	%%% now we can search for those p for which dist(p)=DELTA.
	%%% that means: find zeros of mkx4p(phase,h,p(indy),vp,vs,r,rp)-delta for each piece
	%%% zeros are points at which y=0 is touched or crossed.
	%%% since DELTA might be a vector, we loop over all elements of DELTA.
    deltalen=length(delta);
    for deltacnt=1:deltalen
	   anz=length(starts);
	   for i=1:anz
          %%% search in the following piece
          distpiece=dist(starts(i):ends(i))-delta(deltacnt);
          ppiece=p(starts(i):ends(i));
          anglepiece=angles(starts(i):ends(i));
          piecelen=length(distpiece); % MK19112003, was: ends(i)-starts(i);
          %%% find left an right limits between zeros must be
          %%% (sorry for trickyness)
          indies=find([diff(sign(distpiece)) 0]~=0); % the last index before crossing y=0
          indies=[indies find(abs(distpiece)<epsilondeg)]; % add 0-touches to the list
          indies=unique(indies); % might be non-unique after the previous step!
       
          %%% we now have a list of ray parameters, which are the "last" value before
          %%% a zero. we go through this list and determine a better estimate of the
          %%% location of the zero.
          indyanz=length(indies);
          for indy=1:indyanz
             if piecelen>1 % MK18052005
                 if (indies(indy)+1)<=piecelen
                    %%% probable solution interval is fully covered by ray parm scan
                    [pieceres,a,d]=mkfindzeros([ppiece(indies(indy)) ppiece(indies(indy)+1)],...
                                 [distpiece(indies(indy)) distpiece(indies(indy)+1)]+delta(deltacnt),...
                                 [anglepiece(indies(indy)) anglepiece(indies(indy)+1)],...
                                 phase,h,model,delta(deltacnt),vp,vs,fzopt);
                 else
                    %%% probable solution interval begins at last element of scan
                    %%% searching between last two points may give a good approximation
                    %%% and is probably stable - at least no out of range error occurs.
                    %%% This could probably be improved. MK29.09.2003
                    [pieceres,a,d]=mkfindzeros([ppiece(indies(indy)-1) ppiece(indies(indy))],...
                                 [distpiece(indies(indy)-1) distpiece(indies(indy))]+delta(deltacnt),...
                                 [anglepiece(indies(indy)-1) anglepiece(indies(indy))],...
                                 phase,h,model,delta(deltacnt),vp,vs,fzopt);
                 end; % if (indies(indy)+1)<=piecelen
             else
                 %%% single-sample pieces cannot be analyzed MK18052005
                 pieceres=NaN;
                 d=NaN;
                 a=NaN;
             end; % if piecelen>1
             if abs(d-delta(deltacnt))<epsilondeg
                  pres=[pres pieceres];
                  dres=[dres d];
                  ares=[ares a];
                  deltain=[deltain delta(deltacnt)];
             else % else added MK14102003
                  %%% if no ray parameter was found, return NaN!
                  pres=[pres NaN];
                  dres=[dres NaN];
                  ares=[ares NaN];
                  deltain=[deltain delta(deltacnt)];
             end; % if d-delta<epsilondeg
          end; % for indy
      end; % for i
	end; % for deltacnt
	
	
	
	%%% return result
	p=pres; %(find(~isnan(res)));
	if isempty(p)
	  p=NaN;
	end; % if isempty
	a=ares;
	d=dres; %(find(~isnan(res)));
    
    
    %%% finally, we have to construct the STARTS and ENDS lists for output.
    %%% these range completely different from those determined after
    %%% the continuous scan through all take off angles!
    %%% we therefore have to search for values of a (output takeoff angles)
    %%% which are closest to the start and end values determined above.
    %%% (this was behind the allreals-if before MK18052005)
    [asort,sorter]=sort(a); % sort a by size
    newstarts=[];
    newends=[];
    anz=length(starts);
    for indy=1:anz
        %%% a and angles are assumed to be ordered by size.
        %%% then the first element of indies will point to the smallest element of a
        %%% which is larger than angles(ends(indy)), 
        startindies=find(asort>=angles(starts(indy)));
        if ~isempty(startindies)
           newstarts=[newstarts startindies(1)];
        end; % if ~isempty(startindies)

        %%% a and angles are still assumed to be ordered by size.
        %%% then the last element of indies will point to the largest element of a
        %%% which is smaller than angles(ends(indy)), 
        endindies=find(asort<=angles(ends(indy)));
        if ~isempty(endindies)
           newends=[newends endindies(end)];
        end; % if ~isempty(endindies)
    end; % for anz
    %%% transform index values for ordering by a into index values for ordering by d
    newstarts=sorter(newstarts);
    newends=sorter(newends);

    %%% index lists should be made unique
    newstarts=unique(newstarts);
    newends=unique(newends);

    starts=newstarts;
    ends=newends;
    %%% (this was behind the allreals-if before MK18052005 - end of shifted passage)

else
   %%% no arrivals exist!
   p=NaN; %[];
   a=NaN;
   d=NaN; %[];
   starts=NaN; %MK18052005
   ends=NaN; %MK18052005
   deltain=NaN; %MK18052005
end; % if ~isempty(allreal)




%%% now we're done.
%toctime=toc; disp(['MKFINDP: elapsed time ' mks2hms(toctime)]);


return



