function dmy=mkplotmodel(model)
% MKPLOTND.........plot velocity model read from nd file
%
% call: mkplotnd(model);
%
%      A structure describing the velocity distribution.
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
% result: none
%
% creates several subplots: diagrams of parameters versus depth and cartoons
% that visualize simplified parameter sets in cross sections through the
% planet.
% in these cartoons, the parameter right below a discontinuity is used to
% fill the space between this and the next discontinuity!
% shades show ratio between actual and max value of a parameter.
% white means parameter is zero or undefined.
%
% assumes that the larges depth given in vmod is equal to the plantes radius!
%
% This routine is derived from MKPLOTND, which did the same with a different data
% structure.
%
% Martin Knapmezer, 05.04.2002, 12.04.2002, 05.07.2002


%%% prepare plot window
pos=get(gcf,'Position');
set(gcf,'Position',[pos(1) 0 650 900]);



%%% transform model structure into the old VMOD-array
lanz=length(model.vp); % number of layers
vmod=zeros(lanz,6);
vmod(:,1)=model.z;
vmod(:,2)=model.vp;
vmod(:,3)=model.vs;
vmod(:,4)=model.rho;
vmod(:,5)=model.qp;
vmod(:,6)=model.qs;
dz=[model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz];
[dz,sorter]=sort(dz);
dname=strvcat('Conrad','Moho','Olivine ab','Olivine bg','Olivine Perovskite','CMB','ICB',model.dname); % discontinuity names
dname=dname(sorter,:);
modelname=model.name;
%%% remove non-existing standard discontinuities
indies=find(~isnan(dz));
dz=dz(indies);
dname=dname(indies,:);


%%% add empty layer to define planets surface
vmod=[0 -1 -1 -1 -1 -1; vmod];
dz=[0 dz];
dname=strvcat('surface',dname);


%%%%% some important values
anz=length(dz); % number of discontinuities
lanz=lanz+1; % number of layers - add one for the surface!
rplanet=model.rp; % radius of planet
vsmax=max(vmod(:,3)); % max S velocity
vpmax=max(vmod(:,2)); % max P velocity
rhomax=max(vmod(:,4)); % max density
rhomin=min(vmod(:,4)); % min density
qmax=max(max(vmod(:,[5 6]))); % max Q
qmin=min(min(vmod(:,[5 6]))); % min Q
qsmax=max(vmod(:,6)); % max Qs
qpmax=max(vmod(:,5)); %max Qp





%%%%% subplot: v(z) curves
subplot(4,2,1);
plot(vmod(:,2),vmod(:,1),'k-');
hold on
plot(vmod(:,3),vmod(:,1),'k-');
hold off
%% mark discontinuities
set(gca,'YTick',[dz rplanet]);
set(gca,'YGrid','on');
%set(gca,'XGrid','on');
%% axes and decoration
axis ij
mkaxis('y',[min(vmod(:,1)) max(vmod(:,1))]);
mkaxis('x',[0 15]);
xlabel('velocity [km/s]');
ylabel('depth [km]');
title(['velocity model: ' modelname]);
axis square
box on

%%%%% subplot: density
if rhomax~=0
   subplot(4,2,3);
   plot(vmod(:,4),vmod(:,1),'k-');
   axis ij
   mkaxis('y',[min(vmod(:,1)) max(vmod(:,1))]);
   mkaxis('x',[2 15]);
   xlabel('density \rho [g/ccm]');
   ylabel('depth [km]');
   %% mark discontinuities
   hold on
   for indy=1:anz
      plot([2 15],[1 1]*dz(indy),'k:');
   end; % for indy
   hold off
   %% axes and decoration
   axis square
   box on
end; % if rhomax~=0


%%%%% subplot: Qp, Qs
subplot(4,2,5);
indies=find(vmod(:,5)~=0);
plot(vmod(indies,5),vmod(indies,1),'k-');
hold on
plot(vmod(:,6),vmod(:,1),'k--');
hold off
%% mark discontinuities
hold on
for indy=1:anz
   plot([10 100000],[1 1]*dz(indy),'k:');
end; % for indy
hold off
%% axes and decoration
mkaxis('y',[min(vmod(:,1)) max(vmod(:,1))]);
mkaxis('x',[10 100000]);
xlabel('Qs (--) and Qp (-)');
ylabel('depth [km]');
axis ij
axis square
%axis tight
box on
set(gca,'XScale','log');
set(gca,'XTick',[10 100 1000 10000 100000]);


%%%%% subplot: poisson ratio
subplot(4,2,4);
poisson=zeros(lanz,1);
for i=1:lanz
   poisson(i)=mkpoisson(vmod(i,2),vmod(i,3));
end; % for i
plot(poisson,vmod(:,1),'k-');
hold on;
plot([1 1]*0.25,[0 rplanet],'k:');
hold off;
%% mark discontinuities
hold on
for indy=1:anz
   plot([0 0.5],[1 1]*dz(indy),'k:');
end; % for indy
hold off
%% axes and decoration
mkaxis('y',[min(vmod(:,1)) max(vmod(:,1))]);
set(gca,'XTick',0:0.1:0.5);
set(gca,'YTick',0:2000:6000);
xlabel('poisson ratio');
ylabel('depth [km]');
axis ij
axis square
axis tight
box on



%%%%% subplot: cross section cartoon Vp, Vs, rho
subplot(4,2,2);
n=50; % number of points to use in circle arcs
hold on
for indy=1:anz
   radius=rplanet-dz(indy);
   layerindy=find(vmod(:,1)==dz(indy)); % index of layers around discontinuity
   layerindy=layerindy(2); % second of the two possible solutions
   %% encode values in color
   style='k';
   vscol=1-vmod(layerindy,3)/(vsmax);
   vpcol=1-vmod(layerindy,2)/(vpmax);
   if rhomax~=0
      rhocol=1-vmod(layerindy,4)/(rhomax);
   else   
      rhocol=NaN;
   end; % if isnan
   poiscol=1-poisson(layerindy)/max(poisson);;
   %% plot the circle
   rgb=[1 1 1];
   mkfillcircle(0,0,radius,0,pi/2,n,style,rgb*vscol);
   mkfillcircle(0,0,radius,pi/2,pi,n,style,rgb*vpcol);
   mkfillcircle(0,0,radius,pi,3*pi/2,n,style,rgb*rhocol);
   mkfillcircle(0,0,radius,3*pi/2,2*pi,n,style,rgb*poiscol);
   %% label the different sectors of the cartoon
   text(rplanet*0.8,rplanet*0.9,'Vs');
   text(-rplanet*0.9,rplanet*0.9,'Vp');
   text(-rplanet*0.9,-rplanet*0.9,'\rho');
   text(rplanet*0.7,-rplanet*0.9,'poisson');
end; % for indy
hold off
axis square
axis equal
axis tight
axis off
set(gca,'XTickLabel',' ');
set(gca,'YTickLabel',' ');
set(gca,'TickLength',[0 0]);
%colormap(1-gray(256));
%mkcolorbar('horiz',[0 100],'relative value [%]');


%%%%% subplot: cross section cartoon Qp, Qs
subplot(4,2,6);
n=50; % number of points to use in circle arcs
hold on
for indy=1:anz
   radius=rplanet-dz(indy);
   layerindy=find(vmod(:,1)==dz(indy)); % index of layers around discontinuity
   layerindy=layerindy(2); % second of the two possible solutions
   %% encode Q in color
   style='k';
   qscol=1-log(vmod(layerindy,6))/log(qsmax);
   qpcol=1-log(vmod(layerindy,5))/log(qpmax);
   %% plot the circle
   rgb=[1 1 1];
   mkfillcircle(0,0,radius,0,pi,n,style,rgb*qscol);
   mkfillcircle(0,0,radius,pi,2*pi,n,style,rgb*qpcol);
    %% label the different sectors of the cartoon
   text(-rplanet*0.9,rplanet*0.9,'Qs');
   text(-rplanet*0.9,-rplanet*0.9,'Qp');
   end; % for indy
hold off
axis square
axis equal
axis tight
axis off
set(gca,'XTickLabel',' ');
set(gca,'YTickLabel',' ');
set(gca,'TickLength',[0 0]);


%%%%% subplot: correlation of rho and Vp
subplot(4,1,4);
plot(vmod(:,2),vmod(:,4),'ko');
hold on
plot(vmod(:,2),mkbirch(1000*vmod(:,2))/1000,'k-');
hold off
%% axes and decoration
axis([4 15 2 15]);
xlabel('P velocity [km/s]');
ylabel('density \rho [g/ccm]');
title('\rho as function of Vp (circles), compared to Birch''s Law (line) ');
%axis tight
box on



% return dummy result
dmy=0;