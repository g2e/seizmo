function [p,a,d]=mkfzeronest(phase,delta,h,model,startalpha,startdist);
% mkfzeronest........invoke FZERO to minimize delta(alpha) function
%
% call: [p,a,d]=mkfzeronest(phase,delta,h,model,startalpha,startdist);
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
% This routine uses a modified version of MatLab's fzero() function to find the takeoff angle at
% which a given epicentral distance is reached. Since fzero() accepts only
% functions with a single input parameter, a nested function construction
% is necessary to provide all the input that is needed by MKX4P, especially
% the free variable pair takeoff angle / ray parameter.
%
% The modified version of fzero is mkfzero. The modification is that this
% new version has an emergency break for non-stopping iterations.
%
% the routine uses an iterative approach to guarantee that DELTA is reached
% within a given accuracy (fzero itself only guarantees that the takeoff
% angle leading to that delta is found with given accuracy, and these two
% accuracies are connected via the unknown local function gradient).
%
% Martin KNapmeyer, 15.12.2006, 25.01.2007, 27.02.2007

%%% 25012007 use of mkfzero instead of fzero.
%%%          and do not use TolX values smaller than ABOUTZERO.
%%% 27022007 catch NaN output of MKFZERO


%%% init result
p=[];
a=[];
d=[];


%%% if the alpha interval has a length close to zero, we do not even try to
%%% minimize something - diff(startalpha)==0 means the function is vertical
%%% or the input samples are identical
%%% aboutzero is also used below to avoid ridiculously small TolX values. MK25012007
aboutzero=1e-10;
if diff(startalpha)<aboutzero;
   p=NaN;
   a=NaN;
   d=NaN;
   return;
end; % if diff(startalpha)<eps('double')


%%% epsilon for convergence criterion
%%% the iteration is considered as converged if the epicentral distance
%%% differs from the desired value by less than this
disteps=1e-3; % [deg]


%%% emergency break
%%% This is to prevent infinite iteration MK25012007
maxiter=50;


%%% but minimization procedures usually do not want a tolerance in Y, but
%%% in X, so we have to estimate the X tolerance from the desired Y
%%% tolerance. This is done using the function gradient defined by the
%%% input parameters (very crude... but will be improved iteratively)
alphaeps=disteps/abs(diff(startdist)/diff(startalpha));

%%% initialize bookkeeping
%%% for the iteration of start intervals and tolerances, we have to keep
%%% lists of all fzero results and tolerances, in order to be able to
%%% choose the shortest possible alpha interval for each call of fzero
alphalist=startalpha; % to collect all alpha output of fzero
fvallist=startdist-delta;  % to collect all fval output of fzero
epslist=[]; % to collect all termination tolerances


%%% initialize fzero input
newx=startalpha;


%%% call fzero
%%% we do this in an iterative procedure: since we need an alpha, which
%%% reproduces DELTA within a given tolerance, we repeat fzero with smaller
%%% and smaller input intervals and TolX values until we're good enough
done=0;
evalcnt=0; % objective function evaluations counter
itercnt=0;
while done==0
    
    %%% initialize next fzero call
    %disp(['MKFZERONEST: iteration #' int2str(itercnt+1) ' fzero TolX: ' num2str(alphaeps)]);
    x=newx;
    fzoptions=optimset('Display','off',...
                       'FunValCheck','off',...
                       'MaxIter',maxiter,...
                       'TolX',alphaeps);
    
    %%% execute fzero() with the current settings
    [alpha,fval,exitflag,fzeroinfo]=mkfzero(@mkdisterror,x,fzoptions);
    
    %%% did the search fail?
    %%% MKFZERO may retyurn NaN if it cannot find a solution or if other
    %%% problems occur which require to cancel the search 
    %%% in that case, we have no choice but returning NaN also.MK27022007
    if isnan(alpha)
       p=NaN;
       a=NaN;
       d=NaN;
       return;
    end; % if isnan(alpha)
    
    %%% some bookkeeping
    itercnt=itercnt+1;
    evalcnt=evalcnt+fzeroinfo.funcCount;
    alphalist=[alphalist alpha];
    fvallist=[fvallist fval];
    epslist=[epslist alphaeps];
    
    %%% is fval small enough? if yes, stop, if no: continue
    if abs(fval)<disteps
       %%% we're close enough to zero!
       done=1;
       %disp(['MKFZERONEST: fval=' num2str(fval)]);
    else
       %%% we're not close enough, so we have to do another run of fzero
       %%% with a shorter start interval and a smaller Termination tolerance
       
       %%% determine new x interval
       %%% lower boundary: highest negative function value in fvallist
       %%% uper boundary: smallest positive function value in fvallist
       xnew(1)=max(alphalist(fvallist<0));
       xnew(2)=min(alphalist(fvallist>0));
       
       %%% determine new termination tolerance alphaeps
       alphaeps=min(alphaeps/10,abs(diff(xnew))/10);
       
       %%% alphaeps too small to make sense? MK25012007
       if alphaeps<aboutzero
          done=1;
       end; % if alphaeps<aboutzero
       
    end; % if fval<disteps
    
    %%% iteration exhaust? - just as emergency brake MK25012007
    if itercnt>maxiter
       done=1;
    end; % if itercnt>maxiter
    
end; % while


    %%%%%%%%%%% nested function to be minimized by FZERO
    function y=mkdisterror(x)
    % mkdisterror......compute distance error with respect to target distance
    %
    % call: y=mkdisterror(x);
    %
    %       x: takeoff angle [deg]
    %
    % result: y: difference between desired epicentral distance DELTA as
    %            given in the input of MKFZERONEST and the actual
    %            epicentral distance reached with take off angle ALPHA
    %
    % This routine uses PHASE, DELTA, H, and MODEL from parent function
    % MKFZERONEST as if it were global variables.
    %
    % Martin Knapmeyer, 15.12.2006
    
        %%% evaluate objective function at position X
        rayp=mkangle2rayp(phase,h,x,model); % ray parameter corresponding to XNEU angle
        dist=mkx4p(phase,h,rayp,model,x); % epicentral distance for this angle/rayp
        y=dist-delta; % objective function which is to be zeroed
    
    end % function mkdisterror
    %%%%%%%%%%% end of nested function


% %%% control output
% disp(['MKFZERONEST: function evaluations: ' int2str(evalcnt) ' in ' int2str(itercnt) ' fzero calls.']);

%%% construct output
p=mkangle2rayp(phase,h,alpha,model);
a=alpha;
d=fval+delta;



%%% in nested function constructions, each function has to end with "end"
end