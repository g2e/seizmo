function [p,a,d]=mkfalseposition(phase,delta,h,model,startalpha,startdist);
% mkfalseposition.......determine function root by false position method
%
% call: [p,a,d]=mkfalseposition(phase,delta,h,model,startalpha,startdist);
%
%
%       phase: string containing a single seismic phase name like 'P', 'S',
%              'ScS','PKPdf', etc.
%              Phase names are case sensitive!
%       delta: epicentral distance [deg]
%              might be a vector of epicentral distances
%           h: focal depth [km]
%       model: A structure describing the velocity distribution. Such a
%              structure can be obtined from MKREADND
%  startalpha: two element vector containing the starting values of
%              take off angle alpha [deg]
%              These must be chosen such that the root is between these two
%              values, and startalpha(1)<startalpha(2)
%   startdist: two element vector containing the starting values of
%              epicentral distance, corresponding to STARTALPHA
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
%
% This routine searches the takeoff angle for which phase PHASE, starting 
% at source depth H reaches the epicentral distance DELTA.
% To do so, it searches the root of mkx4p-delta in the takeoff angle
% interval MINALPHA...MAXAPLHA, using the false position method described
% in NUMERICAL RECIPES.
%
% Martin Knapmeyer, 11.12.2006


%%% init result
p=[];
a=[];
d=[];


%%% epsilon for convergence criterion
%%% the iteration is considered as converged if the epicentral distance
%%% differs from the desired value by less than this
disteps=1e-3; % [deg]


%%% initialize iteration loop
x=[startalpha(1); startalpha(2); 0];
y=[startdist(1)-delta; startdist(2)-delta; inf]; % subtract DELTA since we're searching a root!

%%% determine max number of iterations:
%%% the false position method should be faster than bisection, therefore we
%%% use the expected number of iterations that bisection would need as
%%% emergency brake criterion. The problem with this is: we don't care for
%%% how good we know alpha, but we want to know delta with a prescribed
%%% accuracy! So we have to estimate how good we need to know alpha from
%%% that, by estimating the delta-gradient. MK14122006
alphaeps=disteps/abs(diff(startdist)/diff(startalpha));
maxcnt=50; %ceil(log2(abs(diff(startalpha))/alphaeps));
% disp(['MKFALSEPOSITION: alpha epsilon: ' num2str(alphaeps) ', max number of iterations: ' num2str(maxcnt)]);

%%% function evaluation counter
evalcnt=0;

if abs(min(y))<disteps
   %%% the starting interval is already small enough
   %%% unlikely, but possible
   [d,indy]=min(y);
   a=x(indy);
   p=mkangle2rayp(phase,h,a,model);
else
   %%% the starting iterval is not small enough, so start searching
    
    %%% iteration loop
    done=0;
    cnt=0;
    while done==0


        %%% From existing x and y values, guess the root by linear interpolation.
        %%% We consider x as function of y and interpolat x at position y=0.
        y(3)=0;
        x(3)=mk2ptinterp(y(1:2),x(1:2),y(3));


        %%% evaluate objective function at position XNEU
        rayp=mkangle2rayp(phase,h,x(3),model); % ray parameter corresponding to XNEU angle
        [dist,segx,segz]=mkx4p(phase,h,rayp,model,x(3)); % epicentral distance for this angle/rayp
        y(3)=dist-delta; % objective function which is to be minimized
        evalcnt=evalcnt+1; % increment function evaluation counter
        
        
%         %%% control plot
%         figure(2);
%         hold on
%         plot(x(3),max(segz),'r.');
%         hold off
%         figure(1);
%         hold on
%         plot(x(3),dist,'k.');
%         hold off
%         drawnow;


        %%% figure out which of the three points x(1), x(2), x(3) to use: 
        %%% those that bracket the root.
        if abs(y(3))<disteps
           %%% the last point is close enough, stop search
           p=rayp;
           a=x(3);
           d=dist;
           done=1;
           %disp(['MKFALSEPOSITION: iterations carried out: ' int2str(cnt)]);
        else
           %%% between which two points does the sign change occur?
           if sign(y(1))~=sign(y(3))
              %%% points 1 and 3 have different signs, so we throw away #2
              x(2)=x(3);
              y(2)=y(3);
           else
              if sign(y(2))~=sign(y(3))
                 %%% points 2 and 3 have different signs, so we throw away #1
                 x(1)=x(3);
                 y(1)=y(3);
              else
                 %%% all points have the same sign
                 if sum(abs(y))==0
                    %%% same sign for all three is possible if they
                    %%% are all zero! If so, we succeeded!
                    p=rayp;
                    a=x(3);
                    d=dist;
                    done=1;
                 else
                    %%% same sign is also possible if we lost the root.
                    error('MKFALSEPOSITION: root slipped out of bracket!');
                 end;
              end; % if sign(y(2))~=sign(y(3))
           end; % if sign(y(1))~=sign(y(3))
        end; % if y(3)==0



        %%% increment counter
        cnt=cnt+1;

        %%% emergency brake
        if cnt>maxcnt
           if abs(y(3))<disteps*10
               p=rayp;
               a=x(3);
               d=dist;
               %disp(['MKFALSEPOSITION: iteration exhaust, solution quality: ' num2str(abs(y(3)))]);
           else
               p=NaN;
               a=NaN;
               d=NaN;
               %disp(['MKFALSEPOSITION: iteration exhaust, solution quality: NaN']);
           end; % if abs(y(3))<disteps*10
           done=1;
           % disp('MKFALSEPOSITION: iteration exhaust');
        end; % if cnt>maxcnt


    end; % while done==0

end; % if min(y)<disteps


% %%% control output
% disp(['MKFALSEPOSITION: function evaluations: ' int2str(evalcnt)]);

% %%% control plot
% figure(1);
% hold on
% plot(a,d,'gs');
% hold off
% drawnow;


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         HELPER ROUTINES                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
