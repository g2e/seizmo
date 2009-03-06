function [texth,tickh,circh,p4]=mkpolarframe(rp,dz,framewidth,steps,ticklen,dolabel);
% mkpolarframe.....create empty frame for ray diagram, in pola coordinates
%
% call: [texth,tickh,circh]=mkpolarframe;
%       [framewidth,steps,ticklen,dolabel]=mkpolarframe('getdefaults');
%       [texth,tickh,circh]=mkpolarframe(rp,dz);
%       [texth,tickh,circh]=mkpolarframe(rp,dz,framewidth,steps,ticklen,dolabel);
%
%               rp: radius of the diagram [km]
%                   (since this was originally intended to plot seismic
%                   rays within planets)
%
%               dz: discontinuity depths [km]
%                   depths of planetary discontinuities. At the radii
%                   corresponding to these depths, circles will be drawn.
%                   You can also use this to produce radius "ticks"
%
%       framewidth: relative width of the frame that contains the angle
%                   tick marks [percent]
%
%            steps: tickmark intervals [deg]
%                   This is a list that specifies the length of tick marks
%                   Tick marks come in families: e.g. one large one every
%                   90 degrees, smaller ones every 10 degrees, tiny ones
%                   every degree. Or some other combination you like.
%                   The number of entries in STEPS defines how many such
%                   families are generated, and the values in STEPS define
%                   the angular increment of each of these families.
%
%          ticklen: Tick lengths [multiples of FRAMEWIDTH]
%                   For each of the tick families defines in STEPS, this
%                   list defines the length of the ticks in units of the
%                   FRAMEWIDTH
%
%          dolabel: flag that defines if text labels will be printed at the
%                   tick marks or not. if dolabel(i)==1, then the ticks of
%                   the i-th family defined in STEPS will have labels
%                   telling the angle at which the tick is.
%
%          If no input parameters are given, default values will be used.
%
%          In the 'getdefaults'-version, the program returns the default
%          values for FRAMEWIDTH, STEPS, TICKLEN, DOLABEL.
%
% result: texth: handles to the created text objects.
%         tickh: handles to angular tick mark objects
%         circh: handles to circular tick mark objects (the circles
%                defined by DZ, not the framing circles!)
%
%         framewidth -+
%         steps       |
%         ticklen     +--- these correspond to the respective input
%         dolabel    -+    parameters: if you call MKPOLARFRAME with 0 or 2
%                          input arguments, this output tells you what the
%                          used defaults were.
%
% As you see, this is an extremely flexible tool to create polar coordinate
% diagrams.
%
% Martin Knapmeyer, 16.02.2006, 24.05.2006, 17.11.2006


%%% useful constants
radian=pi/180;
fontsize=10;

%%% another constant, I do not remember what it was good for in the
%%% original code... If this is zero, no tick labels are written.
label360=1;

%%% default parameters
dfltframewidth=0.05; % relative width of frame for ticks
dfltsteps=[1 5 10 90]; % tickmark intervals
dfltticklen=[0.3 0.75 1 2]; % tickmark lenghts: tickmarks of interval spteps(1) will habe ticklen(1)
dfltdolabel=[0 0 1 0]; % whether or not to write labels (only in case label360==1)

%%% understand input
nin=nargin;
switch nin
    case {0}
        %%% user provbides no input
        
        %%% radius definitions
        rp=6371; % IASP91 earth radius
        dz=[20 35 410 660 2889 5153.9]; % IASP91 discontinuity depths

        %%% tickmark definitions 
        framewidth=dfltframewidth; 
        steps=dfltsteps;
        ticklen=dfltticklen; 
        dolabel=dfltdolabel;
        
    case {1}
        
        %%% tickmark definitions, to be returned
        texth=dfltframewidth; 
        tickh=dfltsteps;
        circh=dfltticklen; 
        p4=dfltdolabel;
        
        %%% quit routine here!
        return;
        
    case {2}
        %%% user defined radius and discontinuities
        
        %%% tickmark definitions 
        framewidth=0.05; % relative width of frame for ticks
        steps=[1 5 10 90]; % tickmark intervals
        ticklen=[0.3 0.75 1 2]; % tickmark lenghts: tickmarks of interval spteps(1) will habe ticklen(1)
        dolabel=[0 0 1 0]; % whether or not to write labels (only in case label360==1)
       
    case {6}
        
        %%% all parameters given by user, nothing happens here
        
    otherwise
        error(['MKPOLARFRAME: illegal number of input parameters!']);
end; % switch nin

%%% scale tick len by frame width
ticklen=ticklen*framewidth;

%%% maximum ticklength, needed for proper placing outside frame
maxticklen=max(ticklen);
   
   
%%% init new coordinate frame
%cla;

%%% control parameters for plotting circles
n=720; % circle is an (n-1)-sided polygon
phi0=-(2*pi)/n;
phi360=2*pi;

%%% framing circles
mkcircle(0,0,rp,phi0,phi360,n,'k');
hold on;
mkcircle(0,0,rp*(1+framewidth),phi0,phi360,n,'k');
hold off;

%%% tickmarks as angular scale within framing circles
texth=[];
tickh=[];
anz=length(steps); % so many types of labels to draw
hold on
for i=1:anz
    
   %%% construct ticks and set labels
   plotx=[];
   ploty=[];
   for angle=0:steps(i):360;
       
      %%% current tick
      %%% x plot coordinates
      x=rp*cos((180-angle)*radian)*[1 1+ticklen(i)];
      %%% y plot coordinates
      y=rp*sin((180-angle)*radian)*[1 1+ticklen(i)];
      
      %%% add current tick to tick family
      plotx=[plotx NaN x];
      ploty=[ploty NaN y];

      % write labels
      if (dolabel(i)==1)&(angle<360)&(label360==1)
         if angle<100
            strangle=[' ' int2str(mod(fix(180-angle),360))];
         else
            strangle=int2str(mod(fix(180-angle),360));
         end; % if angle<100
         h=text((rp*(1+maxticklen))*cos(angle*radian),...
               (rp*(1+maxticklen))*sin(angle*radian),...
               strangle);
         set(h,'Rotation',angle);
         set(h,'FontSize',fontsize);
         texth=[texth h];
      end; % if dolabel
   end; % for j
   
   %%% plot ticks of this family
   h=plot(plotx,ploty,'k');
   tickh=[tickh h];
   
end; % for i

%%% discontinuities
circh=[];
anz=length(dz); % so many discontinuities
for i=1:anz
   [x,y,h]=mkcircle(0,0,rp-dz(i),phi0,phi360,n,'k');
   circh=[circh h];
   set(h,'UserData',{'depth' num2str(dz(i))});
end; % for i

%%% adjust axes
hold off
axis tight
%axis equal
axis square
%axis auto
axis off
 