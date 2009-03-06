function [rayhandles,texthandles]=mkraydiagram(p1,p2,p3,p4,p5,p6,p7);
% mkraydiagram.......plots nice diagram of ray paths
%
% call: [rayhandles,texthandles]=mkraydiagram(delta,z,type,rp);
%       mkraydiagram(delta,z,type,rp,'addlabel');
%       mkraydiagram(delta,z,type,dz,dname,rp);
%       mkraydiagram(delta,z,type,dz,dname,rp,'addlabel');
%       mkraydiagram(dz,dname,rp);
%       mkraydiagram(dz,dname,rp,'addlabel');
%
%       delta: vector of epicentral distances [deg]
%              describing the "horizontal" component
%              of the ray geometry
%       z: vector of depths [km]
%          descrbing the "vertical" component of ray geometry
%       type: string conntaining a 'P' or 'S' for each ray segment, defining the
%             wave type (as returned by MKX4P)
%       rp: planetary radius [km]
%       dz: list of discontinuity depths
%       dname: list of discontinuity names
%
% result: rayhandles: vector of handles to all generated line objects (rays only)
%         texthandles: vecxtor of handles to all generated text objects (annotations)
%
% This routine creates nice plots of ray paths.
% Three modi are available:
% 1) mkraydiagram(delta,z,type,rp);
%    mkraydiagram(delta,z,type,rp,'addlabel');
%    plot ray into existing diagram - simply plots the ray to the screen,
%    independent of what ios on the screen (YGWYD / you get what you deserve)
% 2) mkraydiagram(dz,dname,rp);
%    mkraydiagram(dz,dname,rp,'addlabel');
%    initialize a new diagram: clears screen & creates polar coordinate system 
%    with angular scale, discontinuities etc.
% 3) mkraydiagram(delta,z,type,dz,dname,rp);
%    mkraydiagram(delta,z,type,dz,dname,rp,'addlabel');
%    initialize new diagram and plot rayn into it
%
% Specifying 'addlabel' as last parameter causes MKRAYDIAGRAM to add angular
% labels around the diameter of the plot.
% In mode 2, labels are set around the whole diameter. In modes 1 and 3, labels
% are set only at points where rays arrive.
%
% rays are plottet piecewise linear between coordinates given in DELTA and Z.
% If you wish to have a better approximiation to the real ray shape, you
% have to use a velocity model with thinner layers.
%
% Martin Knapmeyer, 2002, 16.02.2006
%
% BUGS: - when plotting phases like PS, the last ray segment gets lost, although computed
%         correctly.
%       - when plotting Phases like SKS, the ray color at the CMB is not correct, but
%         timing and ray geometry is correct.

% 16022006: creation of coordinate frame put into separate routine
%           collection of ray segments into single vector


%%% useful constants
radian=pi/180;
fontsize=10;

%%% initialize some arrays
th=[]; % to collect text handles

%%% understand input
nin=nargin;
switch nin
   case {3} %mkraydiagram(dz,dname,rp);
      mode=2;
      dz=p1;
      dname=p2;
      rp=p3;
      type=NaN;
      addlabel=0;
      label360=0;
   case {4} %mkraydiagram(delta,z,type,rp); or mkraydiagram(dz,dname,rp,'addlabel');
      if ischar(p4)
         mode=2;
         dz=p1;
         dname=p2;
         rp=p3;
         addlabel=0;
         label360=1;
         type=NaN;
      else
         mode=1;
         delta=p1;
         z=p2;
         type=lower(p3);
         rp=p4;
         addlabel=0;
         label360=0;
      end; % if ischar
   case {5} %mkraydiagram(delta,z,type,rp,'addlabel');
         mode=1;
         delta=p1;
         z=p2;
         type=lower(p3);
         rp=p4;
         addlabel=1;
         label360=1;
   case {6} %mkraydiagram(delta,z,type,dz,dname,rp);
      mode=3;
      delta=p1;
      z=p2;
      type=lower(p3);
      dz=p4;
      dname=p5;
      rp=p6;
      addlabel=0;
      label360=0
   case {7} %mkraydiagram(delta,z,type,dz,dname,rp,'addlabel');
      mode=3;
      delta=p1;
      z=p2;
      type=lower(p3);
      dz=p4;
      dname=p5;
      rp=p6;
      addlabel=1;
      label360=1;
   otherwise
      error('MKRAYDIAGRAM: use 3, 4, 5, 6, or 7 input parameters!');
 end; % switch nin
 
 if (mode==2)|(mode==3)
     th=mkpolarframe(rp,dz);
 end; % if mode==2 | mode ==3
 
 
 
 %%% plot ray
 anz=length(type); % so many ray segments to plot
 handles=[];
 if ((mode==1)|(mode==3))&(anz>1)
    hold on
    
    %%% empty vectors to colect ray leg plot coordinates
    legx=[];
    legy=[];
    
    %%% loop over all ray segments
    for indy=1:anz;
        
       %%% determine line style from wave type
       switch type(indy)
          case{'p','P'}
             style='b';
          case{'s','S'}
             style='r';
          otherwise
             error(['MKRAYDIAGRAM: type ' type ' unknown.']);
       end; % switch type
       
       %%% plot coordinates of inty-th ray segment are
       x=[(rp-z(indy)).*cos(mod(180-delta(indy),360).*radian)...
          (rp-z(indy+1)).*cos(mod(180-delta(indy+1),360).*radian)];
       y=[(rp-z(indy)).*sin(mod(180-delta(indy),360).*radian)...
          (rp-z(indy+1)).*sin(mod(180-delta(indy+1),360).*radian)];

       %%% add coordinates for this segment to leg list
       legx=[legx x];
       legy=[legy y];
   
             
       %%% if wave type changes, plot current leg and delete coordinates
       %%% from leg point list
       if indy<anz
           if (strcmp(type(indy),type(indy+1))==0) 
               %%% wave type will change from this segment to next one
               %%% OR we at the ray's end
               %%% plot current segment collection with current style and
               %%% delete segement coordinates list
               h=plot(legx,legy,style);
               set(h,'tag',['ray ' type(indy)]);
               handles=[handles h];
               legx=[];
               legy=[];

           else
               %%% wave type will not change
               %%% => continue to collect segments


           end; % (strcmp(type(indy),type(indy-1))==0)
       else
                      
           %%% plot current segment collection with current style and
           %%% delete segement coordinates list
           %disp(['nodes-segments: ' int2str(length(z)-indy)]);
           h=plot(legx,legy,style);
           set(h,'tag',['ray ' type(indy)]);
           handles=[handles h];
           legx=[];
           legy=[];
           
       end; %if indy<anz
       
       
         
%        %%% this plots the ray segment number INDY
%        handles(indy)=plot(x,y,style);
%        set(handles(indy),'tag',['ray ' type(indy)]);
       
    end; % for indy
    
    
    %%% and this adds annotation: epicentral distance of ray end
    if addlabel==1
       [framewidth,steps,ticklen,dolabel]=mkpolarframe('getdefaults');
       ticklen=ticklen*framewidth;
       angle=mod(180-delta(length(delta)),360);
       if ~isnan(angle)
          roundedangle=round(angle*10)/10;
          if angle<100
             strangle=[' ' num2str(180-roundedangle)];
          else
             strangle=num2str(180-roundedangle);
          end; % if angle<100
          h=text((rp*(1+max(ticklen)))*cos(angle*radian),...
                (rp*(1+max(ticklen)))*sin(angle*radian),...
                strangle);
          set(h,'Rotation',angle);
          set(h,'FontSize',fontsize);
          th=[th h];
       end; % if ~isnan(angle)
    end; % if addlabel==true
    
    hold off
    
    
 end; % if mode==1 | mode==3
 
 
 %%% modify axes: labeling must not overwrite axes title
 axis(mkgetaxes(gca));
 
 %return only non-zeros handles, since 0 is the root object...
 rayhandles=handles(find(handles~=0));
 texthandles=th;
 
   