function handles=mkplotttcurves(ttcurves,linespec,xtime);
% mkplotttcurves.........plot travel time curves fro structure
%
% call: mkplottraveltimecurves(ttcurves,linespec);
%       mkplottraveltimecurves(ttcurves,linespec,xtime);
%
%       ttcurves: structure describing travel time curves. Such structures
%                 are produced by MKTTCURVES and consist of the following fields:
%                   .name    : model.name as given by the user
%                   .dangle  : the DANGLE goiven by the user
%                   .h       : focal depth as given by the user
%                   .anz     : number opf phases contained in structure (==length of TTC)
%                   .ttc     : sub-structure of arrays containing the data itself
%                              ttc is an array. each element of ttc describes the travel time curve
%                              of a single phase.
%                   .ttc(i).p: name of the i-th phase
%                   .ttc(i).d: epicentral distances for the i-th travel time curve
%                   .ttc(i).t: travel time for the i-th phase
%       linespec: a LINSPEC string, specifying the plot style.
%          xtime: if present, the x axis will be used as time axis
%
% result: handles: handles to the generated line objects
%
% plots the travel time curves defined by TTCURVES into the current axis.
% SWITCHES HOLD ON AND OFF TO DO SO.
% Distances will be reduced to the interval 0...180deg.
%
% Martin Knapmeyer, 30.07.2002

%%% prepare output
handles=zeros(ttcurves.anz,1);


%%% plot stuff
for indy=1:ttcurves.anz;
   %%% do not attempt if no curve data present!
   if ~isempty(ttcurves.ttc(indy).d)
      %disp(ttcurves.ttc(indy).p);
      %%% reduce distances to 0...180
      ttcurves.ttc(indy).d=mod(ttcurves.ttc(indy).d,360);
      indies=find(ttcurves.ttc(indy).d>180);
      ttcurves.ttc(indy).d(indies)=360-ttcurves.ttc(indy).d(indies);
      %%% and now: plot!
      if nargin==2
          hold on;
          handles(indy)=plot(ttcurves.ttc(indy).d,ttcurves.ttc(indy).t,linespec);
          set(handles(indy),'Tag',ttcurves.ttc(indy).p);
          firstindy=find((~isnan(ttcurves.ttc(indy).d))&(~isnan(ttcurves.ttc(indy).t)));
          if ~isempty(firstindy)
             firstindy=firstindy(1);
             th=text(ttcurves.ttc(indy).d(firstindy),ttcurves.ttc(indy).t(firstindy),...
                     ttcurves.ttc(indy).p,...
                     'Color',get(handles(indy),'Color'),...
                     'VerticalAlignment','bottom',...
                     'Clipping','on');
             textent=get(th,'Extent');
             tpos=get(th,'Position');
             if tpos(1)<90
                tpos(1)=max(tpos(1),0); % positions text at 0deg bundary inside plot
             else
                tpos(1)=min(tpos(1),180-textent(3)); %% positions text at 180 boundary inside plot
             end; % if tpos
             set(th,'position',tpos);
          end; % if ~isempty
          hold off;
      else
          %%% nargin==3: use x-axis as time axis to be compatible with MKPLOTSECTION
          hold on;
          handles(indy)=plot(ttcurves.ttc(indy).t,ttcurves.ttc(indy).d,linespec);
          set(handles(indy),'Tag',ttcurves.ttc(indy).p);
          firstindy=find((~isnan(ttcurves.ttc(indy).d))&(~isnan(ttcurves.ttc(indy).t)));
          if ~isempty(firstindy)
             firstindy=firstindy(1);
             th=text(ttcurves.ttc(indy).t(firstindy),ttcurves.ttc(indy).d(firstindy),...
                     ttcurves.ttc(indy).p,...
                     'Color',get(handles(indy),'Color'),...
                     'VerticalAlignment','bottom',...
                     'Clipping','on');
             textent=get(th,'Extent');
             tpos=get(th,'Position');
             if tpos(2)<90
                tpos(2)=max(tpos(2),0); % positions text at 0deg bundary inside plot
             else
                tpos(2)=min(tpos(2),180-textent(4)); %% positions text at 180 boundary inside plot
             end; % if tpos
             set(th,'position',tpos);
          end; % if ~isempty
          hold off;
      end; % if nargin==2
   else
      disp(['MKPLOTTTCURVES: no data for curve #' int2str(indy) ' (Phase: ' ttcurves.ttc(indy).p ')']);
   end; % if ~isempty(ttcurves.ttc(indy).d) else
end; % for indy


%%% return output
% output is already in the output variable.