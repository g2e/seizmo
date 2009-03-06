function rpd=mkraydepth(p,v,r,h,rp,discons);
% mkraydepth.........compute ray penetration depth IN SPHERICAL EARTH
%
% call: rpd=mkraydepth(p,v,r,h,rp,discons);
%
%           p: ray parameter [s/rad]
%           v: velocities [km/s] within planet at radii given in R
%           r: radii at which velocities are defined [km]
%              R(1) is the largest radius (=the smallest depth)
%              R(1)>=R(2)>=R(3)...
%           h: radius at which the focus is
%           rp: planet radius [km]
%           discons: list of discontinuity radii [km]
%                    This is needed to check if a turning point is at a
%                    discontinutiy, since small numerical errors may lead 
%                    to big problems in those cases.
%
% result: rpd: ray penetration depth: smallest radius (=greatest depth)
%              reached by ray with ray parameter P. At this radius, the
%              ray bends back towards the surface.
%
%              NaN if ray does not bend back.
%
% the velocity is assumed to be piecewise linear: v(r)=v0-g*r
% Piecewise constant models won't have ray vortices inside layers.
%
% possible vortex depths are independent of the takeoff angle of the ray.
% You have to decide if vortex depths are important at all! (for a ray
% which goes upward from the source, vortices below the source are not
% relevant!)
%
% input parameters are in spherical coordinates. The Earth flattening
% transformation is applied internally. This is why the
% planet's radius is needed.
%
% Martin Knapmeyer, 22.04.2002, 12.11.2003-17.11.2003, 24.03.2004,
%                   25.08.2004, 04.10.2006
%                   completely rewritten 06.10.2006
%                   14.11.2006

% 14.11.2006 check for horizontal takeoff
%            source may be at layer bottom


%%% some epsilon
%%% this is in units of slowness in flat earth
%%% if ray parameter corresponds to slowness at source depth within this
%%% unvertainty, it is assumed to be identical. The value is not adjusted carefully,
%%% it just should be small. MK14112006
sloepsilon=1e-6;

%%% this is to determin if a ray turning point depth is at a discontinuity:
%%% yes, if the computed depth deviates by less than this from one of the
%%% listed discontinuities in DISCONS
zprecision=1e-6; % [km]



%%% init result
res=NaN;


%%% v and z have to be row vectors
v=v(:)';
r=r(:)';



%%% Flat Earth Transform for velocity model
[vflat,zflat]=mksfer2flat(v,rp-r,rp);

%%% Flat Earth Transform for ray parameter
pflat=mkpsfer2flat(p,rp);

%%% Flat Earth Transform for source depth
[dmy,hflat]=mksfer2flat(rp-h,rp-h,rp);


%%% convert velocity model into slowness model MK06102006
%%% SLO gives the slowness at which a ray would turn back in depth ZFLAT
%%% The ray turning point is between depth samples i and i+1 if the ray
%%% parameter PFLAT is between slo(i) and slo(i+1)
warning off % to supress warnings that would result from zero velocities
slo=1./vflat;
warning on % but we want other warnings



%%% does the ray parameter mean that we are leaving the source at 90
%%% degrees, i.e. is the vertex at source depth? This is the case if
%%% the ray parameter equals 1/v at source depth in flat earth. MK14112006
%%% And if it is the case, we can return rpd=h and quit here.
sourceslo=mkinterp(zflat,slo,hflat,'all'); % interpolate slowness as source depth
switch length(sourceslo)
    case {0}
         %%% slowness at source depth is undefined!
         %%% We ignore this case here and let the software crash
         %%% elsewhere...
    case {1}
         %%% exactly one slowness value: the question is: has the ray
         %%% parameter the same value? if yes: takeoff==90; rpd=h
         if abs(pflat-sourceslo)<sloepsilon
            rpd=h;
            return;
         end; % if abs(pflat-sourceslo)<sloepsilon
    case {2}
         %%% two slowness values: the source depth is at a discontinuity
         %%% the question is: is the ray parameter within the interval
         %%% defined by these two values? if yes: takeoff==90; rpd=h
         if (pflat>=min(sourceslo))&(pflat<=max(sourceslo))
            rpd=h;
            return;
         end; % 
    otherwise
         %%% I currently have no idea how there could be more than two
         %%% values in a well defined model... but be careful.
         error(['MKRAYDEPTH: unexpected number of slownesses at source depth']);
end; % switch sourceslo




%%% subdivide the slowness model into layers
%%% the i-th layer is between topradius(i) and botradius(i) and has
%%% slowness from topslo(i) to botslo(i).
sampleanz=length(slo); % number of depth samples
topz=zflat(1:sampleanz-1); % upper edges of layers
botz=zflat(2:sampleanz);   % lower edges of layers
topslo=slo(1:sampleanz-1);  % slownesses at upper edges
botslo=slo(2:sampleanz);    % slownesses at lower edges


%%% identify layers which include the current ray parameter
%%% note that slowness decreases with depth if velocity increases with
%%% depth, but slowness increases with depth if velocity decreases with
%%% depth, which may be the case in Low Velocity Zones!
%%% Also note that slowness is constant if velocity is constant.
%%% Therefore we have to check  for slo(top)<=p<=slo(bot) as well as for
%%% slo(top)>=p>=slo(bot)
%%% CANDIDATELAYERS indexes TOPZ and BOTZ and thus defines
%%% layers.
candidatelayers=find(((topslo<=pflat)&(pflat<=botslo))|...
                     ((topslo>=pflat)&(pflat>=botslo)));
                 
                 
if ~isempty(candidatelayers)
   %%% there are candidates for turning depths
                 
                 
    %%% tops and bottoms of candidate layers
    %%% since TOPZ and BOTZ are ordered by depth, TOPCANDIDATE and BOTCANDIDATE
    %%% are also ordered by depth and the first layer in the list is the
    %%% shallowest one.
    topcandidate=topz(candidatelayers); % upper edges of candidate layers
    botcandidate=botz(candidatelayers); % lower edges of candidate layers
    topslocandidate=topslo(candidatelayers); % slowness at candidate layer tops
    botslocandidate=botslo(candidatelayers); % slowness at candidate layer bottoms


%     %%% control plot: slowness model and markers for candidate layers
%     figure(3);
%     clf
%     subplot(1,2,1);
%     plot(slo,zflat,'.-'); % slowness model
%     hold on
%     plot([0 0.4],[1 1]*hflat,'r-'); % source depth
%     plot([1 1]*pflat,[0 45000],'r-'); % ray parameter
%     plot(topslocandidate,topcandidate,'r<'); % top
%     plot(botslocandidate,botcandidate,'r>'); % bottom
%     hold off
%     set(gca,'color',[1 1 1]*0.9);
%     xlabel('Slowness [sec/km]');
%     ylabel('Flat Earth Depth [km]');
%     axis ij
%     subplot(1,2,2);
%     plot(vflat,zflat,'.-'); % velocity model
%     xlabel('Flat Earth Velocity [km/s]');
%     hold on
%     plot([1 max(vflat)],[1 1]*hflat,'r'); % source depth
%     hold off
%     axis ij


    %%% identify layers that are below the source, or contain the source
    deepenough=find(((topcandidate>hflat)&(botcandidate>hflat))|...  % candidate layer below source
                    ((topcandidate<=hflat)&(botcandidate>hflat))|... % source is on layer top or contained in layer
                    ((topcandidate<hflat)&(botcandidate>=hflat))... % source is on layer bottom or contained in layer
                    );    
                
    if ~isempty(deepenough)
        %%% there are layers that are at or below source depth

        %%% remove all others from list
        topcandidate=topcandidate(deepenough);
        botcandidate=botcandidate(deepenough);
        topslocandidate=topslocandidate(deepenough);
        botslocandidate=botslocandidate(deepenough);


%         %%% control plot
%         subplot(1,2,1);
%         hold on
%         plot(topslocandidate,topcandidate,'y<'); % top
%         plot(botslocandidate,botcandidate,'y>'); % bottom
%         hold off


        %%% isolate the shallowest of the candidate layers (which is simply the
        %%% first one in the list, since the list is still ordered)
        topcandidate=topcandidate(1);
        botcandidate=botcandidate(1);
        topslocandidate=topslocandidate(1);
        botslocandidate=botslocandidate(1);


%         %%% control plot
%         subplot(1,2,1);
%         hold on
%         plot(topslocandidate,topcandidate,'g<'); % top
%         plot(botslocandidate,botcandidate,'g>'); % bottom
%         hold off


        %%% interpolate the actual turning point depth from the top and bottom
        %%% slownesses and depths of this layer
        turndepthflat=mk2ptinterp([topslocandidate botslocandidate],...
                              [topcandidate botcandidate],...
                              pflat);

%         %%% control plot
%         subplot(1,2,1);
%         hold on
%         plot(pflat,turndepthflat,'go');
%         hold off
%         subplot(1,2,2);
%         hold on
%         plot([0 max(vflat)],[1 1]*turndepthflat,'g-');
%         hold off


        %%% Inverse Flat Earth Transform for turning point depth
        [dmy,res]=mkflat2sfer(turndepthflat,turndepthflat,rp);
        res=rp-res; % transform returned depth into a radius
        
        
%         %%% control plot
%         subplot(1,2,2);
%         cla;
%         plot(v,rp-r,'.-');
%         hold on
%         plot([0 max(v)],[1 1]*(rp-h),'r-');
%         plot([0 max(v)],[1 1]*(rp-res),'g-');
%         hold off
%         xlabel('Velocity [km/s]');
%         ylabel('Depth [km]');
%         axis ij
%         drawnow
    
    else
        %%% all candidate layers are above the source
        res=NaN;
    end; % if ~isempty(deepenough)

else
   %%% no candidates for tunring depths
   res=NaN;
end; % if ~isempty(candidatelayers)


%%% check if a turning depth is sufficently close to a discontinuity to
%%% assume that it is exactly at a discontinuity MK04122006
%%% RES should contain exactly one value, we do not test this here. But it
%%% might be NaN, in that case, we do not test it.
if ~isnan(res)
   %%% deviation from all of the discontinuity radii
   dr=abs(discons-res);
   %%% check if there is one closer than the epslion
   indy=find(dr<=zprecision);
   if ~isempty(indy)
      %%% RES is less than ZPRECISION from DISCONS(INDY), so we assume it
      %%% is exactly there.
      res=discons(indy);
   end; % if ~isempty(indy)
end; % if ~isnan(res)


%%% return result
rpd=res;





