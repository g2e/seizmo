function [extremalpha,extremp,extremdelta]=mkfindextremedists(alphalimit,deltalimit,phase,h,imodel);
% mkfindextremedists.....identify local extrema of delta(alpha)
%
% call: [extemalpha,extremp,extremdelta]=mkfindextremedists(alphalimit,deltalimit,phase,h,imodel);
%
%       alphalimit: two element vector defining the take off angle
%                   interval to be searched [deg]
%       deltalimit: two element vector defining the epicentral distances
%                   corresponding to ALPHALIMIT[deg]
%
%            phase: string denoting the seismic phase under consideration
%
%                h: focal depth [km]
%
%           imodel: velocioty model as returned by MKIMPROVEMODEL
%
% results: extremalpha: take off angle for which the epicentral distance
%                       becomes locally extremal (min or max!)
%              extremp: ray parameter corresponding to EXTREMALPHA [s/deg]
%          extremdelta: epicentral distance reached with a takeoff angle of
%                       EXTREMALPHA
%
%         all three output variables will be empty if there is no extremum
%         other than the interval ends.
%
% The routine assumes that the delta(alpha) function is well behaved in the
% interval defined by the input, with only one extremum, no poles and no
% holes.
% The routine does not leave the interval defined by the input.
% The search algorithm is the Golden Section Search, but sometimes uses a
% gradient method to determine what type of extremum is searched.
%
% Martin Knapmeyer, 21.11.2006

%%% 19122006 insufficient internal subdivision #1 warning removed (was unnecessary)


%%% init result
extremalpha=[];
extremp=[];
extremdelta=[];

%%% sort input: first takeoff angle must be smaller than second
[alphalimit,sorter]=sort(alphalimit);
deltalimit=deltalimit(sorter);

%%% determine velocity at source depth (needed for angle-ray parm conversion)
%%% if h is a discontinuity depth, there will be two velocity values.
%%% vp(1) (or vs(1)) will be for above discontinuity, vp(2) (or vs(2)) for below.
focus=mkinterpmodel(imodel,h,'simple');
indy=find(focus.z==h);
vp=focus.vp(indy);
vs=focus.vs(indy);


%%% if input points are identical, quit
identepsilon=1e-8; % if input alpha and distance values differ by less than this,
                   % they're assumed to be identical.
if (diff(alphalimit)<identepsilon)&&(diff(deltalimit)<identepsilon)
   extremalpha=alphalimit(1);
   extremdelta=deltalimit(1);
   extremp=mka2p(phase,h,extremalpha,vp,vs,imodel.rp);
   return;
end; % if (alphalimit(1)==alphalimit(2))&(deltalimit(1)==deltalimit(2))


%%% a useful constant
phi=(1+sqrt(5))/2; % golden ratio
goldrat=phi-1; % factor used to subdivide intervals

%%% golden section search terimnation control
gradientstep=0.00001; %0.001; % gradient is measured on such fraction of interval length
goldepsilon=sqrt(eps); % resolution tolerance: minimum bracket width, 
                       % as function of machine epsilon (see NUmerical
                       % Recipes for why it is this value)
%%% maximum number of iterations in golden section part     
%%% this number is based on the fact that the search interval length is
%%% shortened by a factor GOLDRAT in each iteration
goldmaxcnt=ceil((log(alphalimit(end)-log(alphalimit(1)))-log(goldepsilon))/log(phi));
%disp(['MKFINDEXTREMEDISTS: expected number of iterations: ' num2str(goldmaxcnt)]);
goldmaxcnt=goldmaxcnt+1; % to be on the safe side



%%% store all alpha/delta pairs computed during search!
%%% we pay for it, we want it!
alphalist=alphalimit; %%% take off angle list
deltalist=deltalimit; %%% epicentral distance angle list
alphalen=length(alphalist);
%%% generate ray parameter list corresponding to alphalist
rayplist=zeros(size(alphalist));
for alphacnt=1:alphalen
    rayplist(alphacnt)=mka2p(phase,h,alphalist(alphacnt),vp,vs,imodel.rp);
end; % for alphalen



%%% initialize initial bracketing
a(1)=alphalimit(1); % supplied by caller
a(3)=alphalimit(2); % supplied by caller
a(2)=a(1)+goldrat*abs(a(3)-a(1)); % golden ratio subdivision
p(1)=rayplist(1);
p(3)=rayplist(2);
p(2)=mka2p(phase,h,a(2),vp,vs,imodel.rp);

%%% initial evaluation
d(1)=deltalimit(1);
d(3)=deltalimit(2);
d(2)=mkx4p(phase,h,p(2),imodel,a(2));

%%% initial bracketing
%%% we search an angle a(2) with a(1)<a(2)<a(3) for which the epicentral
%%% distance is smaller than both initial values (defines a minimum) or
%%% bigger than both initial values (defines a maximum)
%%% if the a(2) initialized above does not satisfy this, we try other
%%% values until we found one. From the list of angles which we then have
%%% we can choose the smalles interva that bracktets the extremum.
extremumfound=0; % we did not yet find the extremum
done=0;
cnt=0;
while done==0
    
    %%% append new value pair to result lists
    alphalist=[alphalist; a(2)];
    deltalist=[deltalist; d(2)];
    rayplist=[rayplist; p(2)];
    
    %%% control plot
%     subplot(1,2,1);
%     hold on
%     plot(a(2),d(2),'r.');
%     hold off
%     drawnow;
    
    
    %%% check current a(2), d(2) pair
    qualifier=mkqualifyextremum(a,d);
    
    %%% evaluate qualification
    switch qualifier
        case {-1}
             %%% we have a minimum
             searchwhat=qualifier;
             done=1;
        case {+1}
             %%% we have a maximum
             searchwhat=qualifier;
             done=1;
        case {0}
             %%% we need a new middle point
             searchwhat=0;
             
             %%% as first try, we just try subdivision from the other side
             if cnt==0
                a(2)=a(3)-goldrat*abs(a(3)-a(1));
                p(2)=mka2p(phase,h,a(2),vp,vs,imodel.rp);
                d(2)=mkx4p(phase,h,p(2),imodel,a(2));
                %%% add these to the archive lists
                alphalist=[alphalist; a(2)];
                deltalist=[deltalist; d(2)];
                rayplist=[rayplist; p(2)];
             else
                %%% otherwise we try something else... we especially try to be clever.
                %%% We assume that there is only one extremum, either a
                %%% minimum or a maximum. We need to determine which of the
                %%% two we have. To do so, we look at the gradient at the
                %%% ends of the interval. Since an extremum means that the
                %%% sign of the gradient changes, the gradients at the
                %%% interval ends must point into the right direction
                %%% already. While doing this, we think of the two interval
                %%% endpoints as of a (lefth end) and c (right end) MK22112006
                
                %%% first determine the gradients by evaluating a neighbour
                %%% point of each interval end. this must be inside the
                %%% interval and can be re-used later.
                intervallength=sqrt((a(3)-a(1)).^2+(d(3)-d(1)).^2);
                
                rightofa=a(1)+(gradientstep*intervallength); % takeoff angle slighty larger than at a(1)
                leftofc=a(3)-(gradientstep*intervallength); % takeoff angle slightly smaller than at a(3)
 
                prightofa=mka2p(phase,h,rightofa,vp,vs,imodel.rp); % ray parm corresponding to rightofa
                pleftofc=mka2p(phase,h,leftofc,vp,vs,imodel.rp); % ray parm corresponding to leftofc
                
                frightofa=mkx4p(phase,h,prightofa,imodel,rightofa); % function value at rightofa
                fleftofc=mkx4p(phase,h,pleftofc,imodel,leftofc); % function value at leftofc
                
                grada=(frightofa-d(1))/(gradientstep*intervallength); % gradient at a (==a(1))
                gradc=(d(3)-fleftofc)/(gradientstep*intervallength); % gradient at c (==a(3))
                
                %%% add these to the archive lists
                alphalist=[alphalist; rightofa; leftofc];
                deltalist=[deltalist; frightofa; fleftofc];
                rayplist=[rayplist; prightofa; pleftofc];
                
%                 %%% control plot of a and c
%                 subplot(1,2,1);
%                 hold on
%                 text(a(1),d(1),'A','Fontsize',14,'FontWeight','demi');
%                 text(a(3),d(3),'C','Fontsize',14,'FontWeight','demi');
%                 hold off
                
%                 %%% control plot of the gradient defining points
%                 subplot(1,2,1);
%                 hold on
%                 plot([rightofa leftofc],[frightofa fleftofc],'r.');
%                 hold off
                
                %%% now we have to distinguish a lot of cases...
                fdiscr=sign(d(1)-d(3)); % function value discriminant
                switch fdiscr
                    
                    case {-1}
                        %%% f(a)<f(c)
                        if grada>0
                           %%% f'(a)>0
                           if gradc>=0
                              %%% f'(c)>=0
                              %%% point c is the maximum of the interval!
                              searchwhat=+1;
                              extremalpha=a(3);
                              extremp=p(3);
                              extremdist=d(3);
                              done=1;
                              extremumfound=1;
                           else
                              %%% f'(c)<0
                              %%% there is a max between a and c (was insuff interal subdiv #1, revised MK19122006)
                              searchwhat=+1;
                              a(2)=rightofa;
                              p(2)=prightofa;
                              d(2)=frightofa;
                           end; % if gradc>=0
                        else
                           %%% f'(a)<=0
                           if gradc>0
                              %%% f'(c)>0
                              %%% there is a min between a and c
                              searchwhat=-1;
                              a(2)=rightofa;
                              p(2)=prightofa;
                              d(2)=frightofa;
                              done=1;
                           else
                              %%% f'(c)<=0
                               %%% this means: more than one extremum!
                              warning(['MKFINDEXTREMEDISTS: insufficient interval subdivision (#2)!',...
                                       ' (Delta=' num2str(deltalimit(1)) '...' num2str(deltalimit(2)) ')']);
                           end; % if gradc>0
                        end; % if grada
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
                    
                    case {0}
                        %%% f(a)==f(c)
                        switch sign(grada)
                            case {-1}
                                %%% f'(a)<0
                                %%% there is a min between a and c
                                searchwhat=-1;
                                a(2)=rightofa;
                                p(2)=prightofa;
                                d(2)=frightofa;
                                done=1;
                            case {0}
                                %%% f'(a)==0
                                %%% insufficient information! we cannot
                                %%% decide what type of extremum is here,
                                %%% or if there is one at all.
                                %%% (may be it would be useful to start testing random points here)
                                error('MKFINDEXTREMEDISTS: insufficient information!');
                            case {+1}
                                %%% f'(a)>0
                                %%% there is a max between a and c
                                searchwhat=+1;
                                a(2)=rightofa;
                                p(2)=prightofa;
                                d(2)=frightofa;
                                done=1;
                            otherwise
                                error('MKFINDEXTREMEDISTS: unexpected gradient sign (#1)!');
                        end; % switch sign(grada)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
                    
                    case {+1}
                        %%% f(a)>f(c)
                        switch sign(grada)
                            case {-1}
                                %%% f'(a)<0
                                if gradc>=0
                                   %%% f'(c)>=0
                                   %%% there is a minimum between a and c
                                   searchwhat=-1;
                                   a(2)=leftofc;
                                   p(2)=pleftofc;
                                   d(2)=fleftofc;
                                   done=1;
                                else
                                   %%% f'(c)<0
                                   %%% c is the minimum!
                                   searchwhat=-1;
                                   extremalpha=a(3);
                                   extremp=p(3);
                                   extremdist=d(3);
                                   done=1;
                                   extremumfound=1;
                                end; % if gradc>=0
                            case {0}
                                %%% f'(a)==0
                                %%% a is the minimum!
                                searchwhat=-1;
                                extremalpha=a(1);
                                extremp=p(1);
                                extremdist=d(1);
                                done=1;
                                extremumfound=1;
                            case {+1}
                                %%% f'(a)>0
                                %%% there must be a max between a and c
                                searchwhat=+1;
                                a(2)=rightofa;
                                p(2)=prightofa;
                                d(2)=frightofa;
                                done=1;
                            otherwise
                                error('MKFINDEXTREMEDISTS: unexpected gradient sign (#2)!');
                        end; % switch sign(grada)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%   
                    otherwise
                        warning(['MKFINDEXTREMEDISTS: unexpected exception. Gradients: A: ',...
                                num2str(grada) ', C: ' num2str(gradc),...
                                ' (Delta=' num2str(deltalimit(1)) '...' num2str(deltalimit(2)) ')']);
                end; % switch fdiscr
                
                
             end; % if cnt==0
             
             
             
    end; % switch qualifier
    
    %%% increment counter
    cnt=cnt+1;
    
    if cnt>2 % more than 2 loopings do not make any sense here!
       %%% after MACCNT trials, we were not able to identify a min or max.
       %%% this means that one of the initial endpoints is the max and the
       %%% other one is the min. Or both are equal.
       %%% in any case. we do not have to search any further.
       done=1;
       extremumfound=1;
    end; % if cnt>maxcnt
    
end; % while done
%%% we have now bracketed the extremum or shown that there is none.
%%% or we have crashed.

%%% control plot of the initial triple
% subplot(1,2,1);
% %cla;
% hold on
% plot(a,d,'r-');
% hold off
% drawnow;

% %% control plot of a direct angle scan in this interval
% subplot(1,2,1);
% hold on
% for an=min(a):0.1:max(a)
%     dmza=an;
%     dmzp=mka2p(phase,h,dmza,vp,vs,imodel.rp); % corresponding ray parm
%     dmzd=mkx4p(phase,h,dmzp,imodel,dmza); % resulting distance 
%     plot(dmza,dmzd,'y.');    
% end; % for angle
% hold off

%%% sort points by a - in the following it is assumed that they're in
%%% sequence!
[a,sorter]=sort(a);
p=p(sorter);
d=d(sorter);

%%% golden section iteration, if necessary
if extremumfound==0
   %%% here comes the golden section part itself.
   %%% at this point, we have three points a(1), a(2), a(3) and their
   %%% function values d(1), d(2), d(3) such that d(2) is either larger or
   %%% smaller than both d(1) and d(3), depending on the type of extremum
   %%% we are searching. Now we have to construct a new point anew according
   %%% to the rules of the golden section search (see Numerical Recipes)
   %%% and iterate this process.
   done=0;
   cnt=0;
   while done==0
       
       %%% determine interval lengths
       %%% intlen is a two element vector
       %%% intlen(1) is length from a(1) to a(2)
       %%% intlen(2) is len from a(2) to a(3)
       intlen=abs(diff(a));
       
       if sum(abs(intlen))>goldepsilon
           %%% convergence not achieved, continue subdivision
       
           %%% now decide where to place the new point
           if intlen(1)>=intlen(2)
              %%% interval from a(1) to a(2) is the longer one
              %%% new point must be between a(1) and a(2)
              newa=a(1)+goldrat*intlen(1);
           else
              %%% interval from a(2) to a(3) is the longer one
              %%% new point must be between a(2) and a(3)
              newa=a(3)-goldrat*intlen(2);
           end; % if intlen(1)>=intlen(2)
           
           %%% evaluate new point
           newp=mka2p(phase,h,newa,vp,vs,imodel.rp); % corresponding ray parm
           newd=mkx4p(phase,h,newp,imodel,newa); % resulting distance
           
%            %%% control plot of the gradient defining points
%            figure(1);
%            subplot(1,2,1);
%            hold on
%            %plot(a,d,'rd');
%            handle=plot(newa,newd,'m.');
%            %set(handle,'color',rgb(cnt+1,:));
%            hold off
%            drawnow;
%            subplot(2,2,2);
%            hold on
%            handle=plot(cnt,newa,'.');
%            %set(handle,'color',rgb(cnt+1,:));
%            hold off
%            drawnow;
%            subplot(2,2,4);
%            hold on
%            handle=plot(cnt,newd,'.');
%            %set(handle,'color',rgb(cnt+1,:));
%            hold off
% %            drawnow;
% %            plot(a,d,'ro');
% %            hold on
% %            plot(newa,newd,'b+');
% %            hold off
%            drawnow;
           
           %%% add new point to archive of evaluated points
           alphalist=[alphalist; newa];
           rayplist=[rayplist; newp];
           deltalist=[deltalist; newd];
           
           %%% add new point to triplet of old ones, to prepare definition
           %%% of a new triplet
           a=[a newa];
           p=[p newp];
           d=[d newd];
           
           %%% sort points by a
           [a,sorter]=sort(a);
           p=p(sorter);
           d=d(sorter);
           
           %%% identify point in this quadruple, where d is min/max
           switch searchwhat
               case {-1}
                   %%% we're searching a minimum
                   [dmy,extremindy]=min(d);
               case {+1}
                   %%% we're searching a maximum
                   [dmy,extremindy]=max(d);
               otherwise
                   error(['MKFINDEXTREMEDISTS: unexpected SEARCHWHAT: ' num2str(searchwhat)]);
           end; % switch searchwhat
           
           %%% the new triplet consists of the identified min/max point and
           %%% its immediate neighbours, or its two predecessors/successors
           %%% if it is at the end of the interval. MK30112006
           switch extremindy
               case {1}
                   newpoints=[1 2 3];
               case {4}
                   newpoints=[2 3 4];
               otherwise
                   newpoints=(extremindy-1):(extremindy+1);
           end; % switch extremindy
           a=a(newpoints);
           p=p(newpoints);
           d=d(newpoints);
       
       else
           %%% convergence achieved: total length of interval smaller than
           %%% desired resolution
           done=1;
           %disp(['MKFINDEXTREMEDISTS: accuracy ' num2str(sum(abs(intlen))) ' achieved after '  int2str(cnt) ' iterations.']);
       end; % if sum(abs(intlen))>goldepsilon
       
       %%% iteration exhaust
       cnt=cnt+1;
       if cnt>goldmaxcnt
           %%% this is an emergency brake to make the loop stop after some
           %%% considerable time. I hope that the current best of all
           %%% evaluated points is good enough, although it does not meet
           %%% the convergence criterion.
           done=1;
           disp('MKFINDEXTREMEDISTS: iteration exhaust!');
       end; % if cnt>goldmaxcnt
       
   end; % while done
end; % if extremumfound==0


%%% pick the identified extremum from the archive of all avaluations
switch searchwhat
    case {-1}
        [extremdelta,indy]=min(deltalist);
    case {+1}
        [extremdelta,indy]=max(deltalist);
    case {0}
        indy=[];
    otherwise
        error('MKFINDEXTREMEDISTS: unexpected SEARCHWHAT value.');
end; % switch searchwhat
extremalpha=alphalist(indy);
extremp=rayplist(indy);


extremalpha=alphalist;
extremp=rayplist;
extremdelta=deltalist;

        
% %%% control plot
% subplot(1,2,1);
% hold on
% plot(extremalpha,extremdelta,'gs');
% hold off
% drawnow;



return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       HELPER ROUTINES                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=mka2p(phase,h,angle,vp,vs,rp);
% mka2p....wrapper around MKANGLE2RAYP and discontinuity discrimination
%
% call: p=mka2p(phase,h,a,vp,vs,rp);
%
%       phase: string containing seismic phase name, like 'P', 'S', 'ScS', 'PKPdf',...
%              The routione does not test whether your ray partameter makes sense for
%              the phase tyou desire!
%           h: hypocentral depth [km]
%       angle: take off angle at source [deg]
%       vp: p wave velocity at focal depth [km/s]
%       vs: s wave velocity at focal depth [km/s]
%       rp: planetary radius [km]
%
% result: p: ray parameter corresponding to ANGLE [s/deg]
%
% VP and VS may be 2 element vectors if the source is at a discontinuity.
% This routine is intended as wrapper especially for this case: it performs
% the upper side / lower side distinction.
%
% Martin Knapmeyer, 21.11.2006

%%% handle discontinuities, if present
if length(vp)==2
       %%% we have two velocities, which means that the focus is at a
       %%% discontinutiy. Then the velocity to use depends on the takeoff
       %%% angle.
       if angle<=90
          %%% downgoing ray: use lower side velocity
          whichvelocity=2;
       else
          %%% upgoing ray: use upper side velocity
          whichvelocity=1;
       end; % if alphalist(alphacnt)<=90
    else
       %%% there's only one velocity, so the choice is simple.
       whichvelocity=1;
end; % if length(vp)==2

%%% the conversion itself
p=mkangle2rayp(phase,h,angle,...
               vp(whichvelocity),vs(whichvelocity),...
               rp);

           
return;  % mka2p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qualifier=mkqualifyextremum(a,d);
% determines if a minimum or a maximum is given by three points
%
% call: qualifier=mkqualifyextremum(a,d);
%
%       a: three element vector, containing the free variable
%       d: three element vector, d is a function of a
%
% result: qualifier: qualifies if the three function values given by A and
%                    D define a minimum or a maximum.
%                    possible values:
%                    -1: minimum
%                    +1: maximum
%                     0: data not sufficient, try another a(2)
%
% Martin Knapmeyer, 21.11.2006

%%% make sure the ionput is sorted by a
[a,sorter]=sort(a);
d=d(sorter);

if (d(2)<d(1))&&(d(2)<d(3))
    %%% we have a minimum
    qualifier=-1;
else
    if (d(2)>d(1))&&(d(2)>d(3))
        %%% we have a maximum
        qualifier=+1;
    else
        %%% we don't know what we have.
        qualifier=0;
    end; % if (d(2)>d(1))&(d(2)>d(3))
end; % if (d(2)<d(1))&(d(2)<d(3))

return; % mkqualifyextremum(a,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
