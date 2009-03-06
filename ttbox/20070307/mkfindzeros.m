function [p,a,d]=mkfindzeros(prange,drange,arange,phase,h,model,delta,vp,vs,options);
% mkfindzeros.........finds zeros of mkx4p(phase,h,p(indy),model)-delta
%
% call: [p,a,d]=mkfindzeros(prange,drange,arange,phase,h,model,delta,vp,vs);
%       [p,a,d]=mkfindzeros(prange,drange,arange,phase,h,model,delta,vp,vs,options);
%
%         prange: Two element vector defining the
%                 range of ray parameters in which to search
%                 It is assumed that exactly one zero lies between prange(1)
%                 and prange(2)
%         drange: epicentral distance (delta) values corresponding to prange values.
%                 used to make a first guess for the solution by linear interpolation
%         arange: take off angle range corresponding to PRANGE
%                 necessary to distinguish between upgoing and downgoing rays
%         phase:  seismic phase name sting, like 'P', 'S', 'ScS' etc.
%                 Phase names are case sensitive!
%             h:  focal depth [km]
%         model:  A velocity model structure, as returned by MKREADND
%         delta:  Epicentral distance for which the takeoff angle and ray
%                 parameter are searched. [deg]
%         vp: p velocity at source depth (may be two elements vector at discontinuities)
%         vs: s velocioty at source depth (may be two elements vector at discontinuities)
%
%         options: optional control parameters affecting the loop control
%                  This is a structure containing the following fields:
%
%                  .initwidth: initial step width in fractions of the
%                              initial search interval width
%                              default is 1.5%, but this shoudl be chosen
%                              depending on the quality of the initial
%                              guess.
%
%                  .wdtdivisor: factor by which step width is divided
%                               when search direction is flipped
%                               default is 1/(0.5*(sqrt(5)-1))=1.618...
%
%                  .maxcnt:    maximum number of iterations
%
%                  .epsilon:   the target function is considered to be zero
%                              when its abs() is smaller than epsilon
%
% result: p: vector containing the ray parameter within PRANGE for which
%            mkx4p returns the DELTA desired by MKFINDP.
%         a: take off angle of the PHASE at the source
%            necessary to distinguish between upgoing and downgoing rays       
%         d: The DELTA that MKX4P returns for P.
%
%         p,a will be empty if no solution exists.
%
% the function actually minimizes abs(mkx4p()-delta) in order to find points
% where y=0 is not crossed but touched.
%
% Algorithm: starts at ARANGE(1) and increases a by some da. flips direction
% as soon as the target function starts running away from zero, and decreases da by some factor.
% This is an iterative process which should converge to the zeros of the target function.
% If there is more than one zero within the ARANGE-interval, it depends on their distribution
% which one is returned. So make sure that your intervals contain only one solution.
% To get a first guess for the soluition, linear interpolation between the
% borders of the interval is used.
%
% Martin Knapmeyer, 07.05.2002, 22.09.2005



%%% understand input
nin=nargin;
switch nin
    case {9}
        initwidth=0.015; 
        wdtdivisor=1/(0.5*(sqrt(5)-1)); 
        epsilon=0.001; % "zero" is when abs()<epsilon
        maxcnt=50; % no more than MAXCNT iterations
    case {10}
        initwidth=options.initwidth; 
        wdtdivisor=options.wdtdivisor; 
        epsilon=options.epsilon;
        maxcnt=options.maxcnt;
    otherwise
        error('MKFINDZEROS: illegal number of input parameters!');
end; % switch nin


%%% initialize search
direct=1; % +1: increase p, -1: decrease p
da=(arange(2)-arange(1))*initwidth; % initial stepwidth for increasing a. may be small since the first guess will be good.
cnt=0; % counter for  iterations
evalcnt=0; % counter for function evaluations
currenta=mk2ptinterp(drange,arange,delta); %currenta=interp1(drange,arange,delta);
if (currenta>=arange(1))&&(currenta<=arange(2)) %if ~isnan(currenta)
   %%% interpolation found approx solution
   if currenta>=90
         whichangle=1;
      else
         whichangle=min(2,length(vp));
      end; % if angle
   currentp=mkangle2rayp(phase,h,currenta,vp(whichangle),vs(whichangle),model.rp);
   %disp(['MKFINDZEROS: prange:' num2str(prange) ' drange: ' num2str(drange) ' 1st guess: ' num2str(currentp)]);
   badness=abs(mkx4p(phase,h,currentp,model,currenta)-delta);
   evalcnt=evalcnt+1;
   oldbad=badness;
   done=0;
%clf;
   while done==0
      if cnt>=maxcnt
         done=1; % give up
         %disp(['MKFINDZEROS: give up for delta=' num2str(delta) '. Current ray parm: ' num2str(currentp)]);
      else
         %%% compute new function value
         oldbad=badness;
         d=mkx4p(phase,h,currentp,model,currenta);
         badness=abs(d-delta);
         
% if cnt>30
%     %figure(1); hold on; plot(currenta,badness,'k.'); text(currenta,badness,[' ' int2str(cnt)]); hold off;
%     figure(1); 
%     hold on;
%        plot(arange,[1 1]*delta);
%        plot(currenta,d,'k.');
%        text(currenta,d,[' ' int2str(cnt)]);
%     hold off;
% end; % if cnt>10
     
         evalcnt=evalcnt+1;
         
         if badness<epsilon
            done=1; % we're close enough!
         else
         
            %%% flip direction if necessary
            %%% interval boundary crossed?
            if (currenta<arange(1))|(currenta>arange(2))
               %disp('MKFINDZEROS: flip due to boundary crossing');
               direct=-direct;
               
            else
                %%% badness got worse?
                if (oldbad<badness)
                   %disp('MKFINDZEROS: flip due to badness increase!');
                   direct=-direct;
                   da=da/wdtdivisor;
                end; % if
            end; % if (currenta<arange(1))|(currenta>arange(2))
            
            %%% new CURRENTA
            currenta=currenta+direct*da;
            
            %%% new CURRENTP
            if currenta>=90
               whichangle=1;
            else
               whichangle=min(2,length(vp));
            end; % if angle
            currentp=mkangle2rayp(phase,h,currenta,vp(whichangle),vs(whichangle),model.rp);
           
            %%% increase iteration count
            cnt=cnt+1; % next iteration...
         end; % if badness else
      end; % if cnt>=maxcnt
   end; % while done
% disp(['MKFINDZEROS: ' int2str(evalcnt) ' evaluations in ' int2str(cnt) ' iterations.']);
else
   %%% interpolation did not find approx. solution - return empty
   currentp=[];
   currenta=[];
   d=[];
end; % if ~isempty(currenta)

% disp(['MKFINDZEROS: iteration count upon exit: ' int2str(cnt)]);
% if cnt==50
%    disp(cnt);
% end; % if cnt
%disp(['MKFINDZEROS: results: p=' num2str(currentp) ', a=' num2str(currenta) ', d=' num2str(d)]);

% figure(2);
% hold on
% plot(delta,cnt,'gx');
% hold off
% drawnow;
% figure(1);


%%% return result
p=currentp;
a=currenta;
% d=d; % this is already done

return;