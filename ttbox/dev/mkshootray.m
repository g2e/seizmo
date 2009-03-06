function [p,a,d,deltain]=mkshootray(phase,deltalist,h,model,pscan);
% mkshootray.........finds ray parameter corrsponding to an epicentral distance IN SPHERICAL EARTH
%
% call: [p,a,d,deltain,starts,ends]=mkshootray(phase,deltalist,h,model);
%       [p,a,d,deltain,starts,ends]=mkshootray(phase,deltalist,h,model,pscan);
%
%       phase: string containing a single seismic phase name like 'P', 'S',
%              'ScS','PKPdf', etc.
%              Phase names are case sensitive!
%   deltalist: epicentral distance [deg]
%              might be a single distance or a vector of epicentral distances
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
%          pscan: output of mkpsampler(phase,h,model);
%                 If these parameters are not given, MKPSAMPLER will be called
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
%                  p(i), a(i), d(i) are the results for input value
%                  deltain(i). This list can be used to check the accuracy
%                  of the indentified d values.
%
% This routine replaces the older MKFINDP routine. MKSHOOTRAY uses the more
% advanced model analysis implemented in MKPSAMPLER, rather than doing a
% brute force attack like MKSCANP did. The new method is faster, more
% robust in findig solutions near critical distances and hits the aimed
% distance with better precision.
% The calling sequence is the same as for MKFINDP. It shoul dbe easy to
% update all programs that do an aimed shooting.
%
% Due to limiiations of zero finding algorithms, MKSHOOTRAY may have
% problems to find rays that go to distances at which delta(alpha) has a
% local minimum (because mkx4p-delta does not change sign there).
%
% Martin Knapmeyer, 11.12.2006, 15.12.2006, 19.12.2006

%%% 15122006 use of MKFZERONEST: neste function construct to call FZERO
%%% 19122006 handle lists of target distances instead of only scalars

%%% initialize result
pres=[];
ares=[];
dres=[];
p=[];
a=[];
d=[];
deltain=[];


%%% some constants
radian=pi/180;
epsilondeg=0.01; % epsilon used in delta-search [deg]


%%% call MKPSAMPLER if necessary. The caller may provide a PSCAN structure,
%%% but if not, we have to build one.
if nargin==4
    p=NaN;
    pscan=mkpsampler(phase,h,model);
end; % if nargin==4


% %%% control plot
% figure(1);
% clf;


%%% loop over all entries in DELTA to handle multiple target distances
deltaanz=length(deltalist);
for deltacnt=1:deltaanz;
    
    %%% current delta to be avaluated
    delta=deltalist(deltacnt);


    %%% loop over the take off angle intervals defined in the PSCAN structure
    %%% by comparing pscan.dist with the requested delta it is easy to select
    %%% those intervals that contain a solution
    %%% an angle interval is defined by the i-th and the (i+1)-th takeoff
    %%% angle.
    angleanz=length(pscan.angles); % number of take off angle samples in PSCAN
    for anglecnt=1:angleanz-1

    %     %%% control plot
    %     figure(1);
    %     plot(pscan.angles,pscan.dist,'b.:');
    %     hold on
    %     plot([0 max(pscan.angles)],[1 1]*delta,'g');
    %     plot(pscan.angles([anglecnt anglecnt+1]),pscan.dist([anglecnt anglecnt+1]),'ro');
    %     hold off
    %     grid on
    %     %axis([14.00769106383256  14.00771690822060  89.60463534663074  89.85063094231177]);
    %     drawnow;

        if pscan.dist(anglecnt)==delta
           %%% interval start sample is exactly at the desired distance, append
           %%% it to the list of solutions (unlikely, but possible)
           ares=[ares; pscan.angles(anglecnt)];
           pres=[pres; pscan.p(anglecnt)];
           dres=[dres; pscan.dist(anglecnt)];
           deltain=[deltain; delta];
        else
            %%% may be the interval end is a solution?
            if pscan.dist(anglecnt+1)==delta
                %%% interval end sample is exactly at the desired distance,
                %%% append it to the list of solutions (unlikely, but possible)
                ares=[ares; pscan.angles(anglecnt+1)];
                pres=[pres; pscan.p(anglecnt+1)];
                dres=[dres; pscan.dist(anglecnt+1)];
                deltain=[deltain; delta];
            else
                %%% the ends of the interval are both no solutions. So we have
                %%% to search for solutions in between. If there are any.
                %%% here we assume that the delta(alpha) function is monotonous
                %%% within the given angle interval!
                if sign(pscan.dist(anglecnt)-delta)~=sign(pscan.dist(anglecnt+1)-delta)
                   %%% there must be a solution in the interval, because the
                   %%% two differences in the if clause have different sign
                   %%% now we search a root of the function dist-delta

                   minimizer='mkfzeronest'; % recommended minimizer is MKFZERONEST MK15122006
                   switch minimizer
                       case {'mkfindzeros'}

                           %%% use MKFINDZEROS
                           %%% THIS SECTION IS DEPRECATED BECAUSE OF BEING UNSTABLE
                            fzopt.initwidth=0.25;
                            fzopt.wdtdivisor=1/(0.5*(sqrt(5)-1));
                            fzopt.epsilon=0.001;
                            fzopt.maxcnt=50;
                            prange=[pscan.p(anglecnt) pscan.p(anglecnt+1)];
                            arange=[pscan.angles(anglecnt) pscan.angles(anglecnt+1)];
                            drange=[pscan.dist(anglecnt) pscan.dist(anglecnt+1)];
                            [p,a,d]=mkfindzeros(prange,drange,arange,...
                                                phase,h,model,delta,...
                                                pscan.vp,pscan.vs,fzopt);
                            ares=[ares; a];
                            pres=[pres; p];
                            dres=[dres; d];
                            deltain=[deltain; delta];
                           %%% end of MKFINDZEROS use


                       case {'mkfzeronest'}

                            %%% use MatLab's FZERO in a nested function construction
                            %%% MK15122006
                            [p,a,d]=mkfzeronest(phase,delta,h,model,...
                                               [pscan.angles(anglecnt) pscan.angles(anglecnt+1)],...
                                               [pscan.dist(anglecnt) pscan.dist(anglecnt+1)]);
                            ares=[ares; a];
                            pres=[pres; p];
                            dres=[dres; d];
                            deltain=[deltain; delta];
                            %%% end of FZERO usage

                       case {'mkfalseposition'}

                            %%% use MKFALSEPOSITION
                            [p,a,d]=mkfalseposition(phase,delta,h,model,...
                                                   [pscan.angles(anglecnt) pscan.angles(anglecnt+1)],...
                                                   [pscan.dist(anglecnt) pscan.dist(anglecnt+1)]);
                            ares=[ares; a];
                            pres=[pres; p];
                            dres=[dres; d];
                            deltain=[deltain; delta];
                            %%% end of MKFALSEPOSITION use

                       otherwise
                           error(['MKSHOOTRAY: unknown function minimizer ' upper(minimizer)]);
                   end; % switch minimizer



                else
                   %%% there is no solution in the interval.
                   %%% in this case, we don't do anything.
                end; % 
            end; % 
        end; % if pscan.dist(anglecnt)==delta

    end; % for anglecnt

    % %%% control output
    % disp('MKSHOOTRAY: achieved accuracy:');
    % disp(num2str(abs(dres-delta)));


    % %%% control plot of results
    % pscanclassic=mkscanp(phase,h,model,0.1);
    % figure(1);
    % hold on
    % plot(pscanclassic.angles,pscanclassic.dist,'k');
    % plot([0 90],[1 1]*delta,'r');
    % plot(ares,dres,'gs');
    % hold off

end; % for delta


%%% return results
if isempty(pres)
    p=NaN;
else
    p=pres;
end; % if isempty(pres)
a=ares;
d=dres;