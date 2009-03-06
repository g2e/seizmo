function handle=mkttdiagram(p1,p2,p3,p4);
% mkttdiagram........draw travel time curves with min and deg scales
%
% call: handle=mkttdiagram(trange);
%       handle=mkttdiagram(delta,tt,style);
%       handle=mkttdiagram(trange,delta,tt,style);
%
%           trange: time range 
%                   maximum traveltime to be displayed in minutes
%                   the diagrams y-axis always begins at t=0
%           delta: vector of epicentral distances [deg]
%           tt: vector of traveltimes [sec]
%               TIMES ARE IN SECONDS as returned by MKTIM4P!
%           style: line style specifier as for PLOT (a linespec)
%
% result: handle: vector of handles to the crated line objects
%                 (only for the travel time curves, not for the 
%                 decoration!)
%
% The routine does not domuch more than call plot(delta,tt,style);
% Most of the work done is for creating a min-scale on the y-axis.
%
% While plotting, HOLD is ON.
% When called with TRANGE argument only, routine prepares an empty
% diagram with time axis from trange(1) to trange(2).
% When called with delta and tt and style, the travel time curve
% is plotted into the current graph, independent of what this graph is.
% When called with all four arguments, a diagram is created and the
% traveltimes plotted into that diagram.
%
% creation of diagrams clears the current axes (cla).
%
% Martin Knapmeyer, 15.05.2002

%%% result
dmy=0;


%%% patch color for axes
patchcol='w';

%%% understand input
nin=nargin;
switch nin
   case {1}
      mode=1;
      trange=sort([0 p1])*60; %%% transform trange from minutes to seconds
   case {3}
      mode=2;
      delta=p1;
      tt=p2;
      style=p3;
   case {4}
      mode=3;
      trange=sort([0 p1])*60; %%% transform trange from minutes to seconds
      delta=p2;
      tt=p3;
      style=p4;
   otherwise
      error('MKTTDIAGRAM: wrong number of input arguments!');
end; % switch nin
   

%%% plot something
switch mode  
   case {1}
      %%% prepare diagram
      %%% MatLab's axes will be scaled in seconds, the minute axis in the plot therefore
      %%% needs some calculations
      %%% general tickmark definitions
      framewidth=5;
      %%% tickmark definitions for minute axis [given in seconds]
      msteps=[1 5 10 60]; % [sec]
      mticklen=[0.3 0.75 1 1]*framewidth;
      mdolabel=[0 0 1 0];
      %%% tickmarks definitions for seconds axis
      ssteps=[10 50 100 500]; % [sec]
      sticklen=[0.3 0.75 1 1]*framewidth;
      sdolabel=[0 0 0 1];
      %%% tickmarks definitions for distance axis
      dsteps=[1 5 10 20]; % [deg]
      dticklen=[0.3 0.75 1 1]*framewidth;
      ddolabel=[0 0 0 1];
      %disp('MKTTDIAGRAM: prepare diagram');
      %%% init axes object
      if ~ishold
         cla;
      end; % if ~ishold
      minx=-framewidth;
      maxx=180+framewidth;
      mint=trange(1);
      maxt=trange(2);
      axis([minx maxx mint-(maxt-mint)*framewidth/1500 maxt+(maxt-mint)*framewidth/1500]);
      %%% draw minutes axis
      patch([minx 0 0 minx],[mint mint maxt maxt],patchcol);
      anz=length(msteps);
      hold on
      for i=1:anz % so many types of labels
         for ttick=mint:(msteps(i)*60):maxt
            x=mticklen(i)*[-1 0];
            y=[1 1]*ttick;
            plot(x,y,'k');
            if mdolabel(i)==1
               lblstr=[num2str(ttick/60) ' '];
               text(x(1),y(1),lblstr,'HorizontalAlignment','right','FontSize',8);
            end; % if mdolabel(i)==1
         end; %for ttick
      end; % for i
      hold off
      %%% draw seconds axis
      patch([maxx 180 180 maxx],[mint mint maxt maxt],patchcol);
      anz=length(ssteps);
      hold on
      for i=1:anz % so many types of labels
         for ttick=mint:ssteps(i):maxt
            x=180+sticklen(i)*[0 1];
            y=[1 1]*ttick;
            plot(x,y,'k');
            if sdolabel(i)==1
               lblstr=[' ' num2str(ttick)];
               text(x(2),y(2),lblstr,'HorizontalAlignment','left','FontSize',8);
            end; % if mdolabel(i)==1
         end; %for ttick
      end; % for i
      hold off
      %%% draw distance axes
      boxhgt=(maxt-mint)*framewidth/200;
      patch([0 180 180 0],[mint mint mint-boxhgt mint-boxhgt],patchcol);
      patch([0 180 180 0],[maxt maxt maxt+boxhgt maxt+boxhgt],patchcol);
      anz=length(dsteps);
      hold on
      for i=1:anz % so many types of labels
         for dtick=(minx+framewidth):dsteps(i):(maxx-framewidth)
            x=[1 1]*dtick;
            y=mint+[0 -1]*dticklen(i)*boxhgt/framewidth;
            plot(x,y,'k');
            %%% write labels to lower axis
            if ddolabel(i)==1
               lblstr=[num2str(dtick)];
               y=mint+[0 -2]*dticklen(i)*boxhgt/framewidth;
               text(x(2),y(2),lblstr,'HorizontalAlignment','center','FontSize',8);
            end; % if mdolabel(i)==1
            %%% plot upper axis
            x=[1 1]*dtick;
            y=maxt+[0 1]*dticklen(i)*boxhgt/framewidth;
            plot(x,y,'k');
         end; %for ttick
      end; % for i
      hold off
      %%% additional settings
      
      ax=mkgetaxes(gca);
      axis(ax);
      axis off
      
      text(90,...
           mint-3*max(dticklen)*boxhgt/framewidth,...
           'Epicentral Distance [deg]',...
           'HorizontalAlignment','center');
      
      text(-framewidth*4,...
           mint+(maxt-mint)/2,...
           'Travel Time [min, sec]',...
           'Rotation',90,...
           'HorizontalAlignment','center');
   case {2}
      %%% draw traveltime curve
      %disp('MKTTDIAGRAM: plotting times');
      %%% transform delta from infinite angle range to 0...180
      %%% with reflection at 180 (this is not simply a mod(delta,180) !
      delta=mod(delta,360);
      indies=find(delta>180);
      delta(indies)=360-delta(indies);
      %%% and now the plotting
      hold on
      handle=plot(delta,tt,style);
      hold off
   case {3}
      %%% draw diagram and curve
      %%% generate diagram by recursive call
      mkttdiagram(trange);
      %%% draw curve by recursive call
      handle=mkttdiagram(delta,tt,style);
end; % switch mode

   
%%% return result
%  handle=handle already   
   