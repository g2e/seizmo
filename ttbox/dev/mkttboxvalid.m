function [h1,h2r,h2tt,h3]=mkttboxvalid(phase,h,dz,testfile,reference,mode,linespec,takeabs,radius);
% mkttboxvalid.......validation of TTBOX calulation against reference data
%
% call: handle=mkttboxvalid(phase,h,dz,testfile,reference,linespec);
%       handle=mkttboxvalid(phase,h,dz,testfile,reference,linespec,takeabs,radius);
%
%       phase: name of phase for which the validation is carried out
%              (although TTBOX does generally not distinguish PKPab and PKPbc,
%               this may be 'PKPab' or 'PKPbc'. But use PKIKP instead of PKPdf.
%               PKPab and PKPbc validations work only for the IASP91 model.)
%           h: focal depth [km]
%          dz: depth sampling interval [km]
%              the model used in computation will be sampled in stepd od DZ km
%       testfile: full path to a CLR or ND file containing the velocity model.
%                 If this is a CLT file, it will be interpolated according to DZ.
%                 If this is a ND file, it will be taken as is.
%                 Radius modification is currently not possible with .nd files.
%       reference: full path to a reference data file
%                This has to be an ASCII file with three columns:
%                1) Epicentral distance [deg]
%                   (it is implicitly assumed that only integer distances are used)
%                2) minutes of travel time
%                3) seconds of travel time
%                (This format is owed to the column format of the IASP91 booklet by Kennett)
%                PLEASE MAKE SURE THAT THESE DATA ARE VALID FOR THE MODEL AND FOCAL DEPTH
%                THAT YOU SPECIFIED IN CLRFILE AND H !! (otherwise the result will be nonsense)
%       mode: model sampling mode of MKCLR2MODEL:
%             'spherical' depth sampling will be equidistant in sphercial coord.
%             'flat' depth sampling will be equidistant in flat coord., which gives
%                    better results but takes longer.
%       linespec: line specifier as used in the plot command.
%       takeabs: string specifying whether or not to plot a signed
%                difference or an absolute difference.
%                possible values: 'abs': plot abs difference
%                                 'sign': plot signed difference
%                                 'rayp': plot signed difference as function of 
%                                         ray parameter (ray parameter will be
%                                         the one found by TTBOX for the distance
%                                         in the reference file)
%                                 'raydepth' plot signed difference as function of
%                                         rays's deepest point (does not work with PcP)
%                                         ray turning points are computed for Vp.
%                DEFAULT: 'abs'
%       radius: this parameter allows to specify a planet radius other than
%               given in the CLR file. RADIUS is a radius in km.
%               Radius modification is currently not possible with nd files.
%               Default: radius as given in CLR file
%
% result: handles to the created line objects:
%           h1: handle to line in figure 1
%          h2r: handle to reference data in figure 2
%         h2tt: handle to TTBOX data in figure 2
%           h3: handle to line in figure 3
%
% An example for reference data lines (valid for IASP91 P phase):
%
% 23 5 6.34
% 24 5 16.31
% 25 5 25.43
% 26 5 34.51
% 27 5 43.54
% 28 5 52.50
%
% These few lines define the P travel time art distances 23-28 deg to be between 5min6.34s
% and 5min52.50sec.
%
% TTBOX will compute travel times for the epicentral distances given in the reference file.
% In case of triplications or other ambuguities, the earliest arrival time will be used.
% The the absolute difference between computed travel time and the travel time given in the
% reference will be plotted.
% Only epicentral distances beyond 10deg are evaluated.
%
% The routine produces three plots:
% figure 1: abs value of travel time deviation, plotted using LINESPEC
% figure 2: travel time curves as computed and as given in reference file. The reference
%           data is always plotted in grey. The TTBOX result is plotted using LINESPEC,
%           but additionally with markers at the positions of the computed points.
% figure 3: round off error of ray parameter estimation, in terms of epicentral distance
%           this is the difference between the DELTAOUT and DELTAIN output of MKFINDP
%           plotted using stem() and LINESPEC on semilog axes.
% The Tag-property of the line objects is set to ['dz=' num2str(dz)] to make later
% identification easier.
%
% All plots are done with HOLD ON. You probably wish to combine the output of several tests
% into one plot.
%
% Martin Knapmeyer, 19.11.2003, 16.03.2004, 29.03.2004, 30.03.2004,
%                   04.12.2006, 19.12.2006

%%% takeabs and radius parameter introduced, MK 15.03., 16.03.2004
%%% distinction of PKPab, PKPbc introduced
%%% evaluation of ND files 30.03.2004
%%% use of MKDETECTDISCON
%%% 19.12.2006 use of MKSHOOTRAY

tic;

%%% init results MK19122006
h1=NaN;
h2r=NaN;
h2tt=NaN;
h3=NaN;

%%% understand input 
if nargin==7
   takeabs='abs';
else
   takeabs=takeabs;
end; % if nargin==7

switch phase
   case {'PKPab'}
      fullphase=phase;
      phase='PKP';
   case {'PKPbc'}
      fullphase=phase;
      phase='PKP';
   otherwise
      fullphase=phase;
end; % switch phase

%%% load reference file
file=load(reference);

%%% decompose reference into distance, min, sec columns
delta=file(:,1);
minutes=file(:,2);
seconds=file(:,3);

%%% compute travel time in sec
tt=minutes*60+seconds;
%%% now plot(delta,tt); would plot the reference travel time curve.


%%% load velocity model
[dmy1,dmy2,extension]=fileparts(testfile);
switch extension
   case {'.clr'}
      clr=mkreadclr(testfile,'silent');
      %%% modify radius?
      if nargin>7
         clr.rp=radius;
      end; % if nargin==7
      %%% discretize continuous model at desired sampling interval
      model=mkclr2model(clr,dz,mode);
   case {'.nd'}
      %%% note that the dz parameter is ignored here!
      model=mkreadnd(testfile,'silent');
   otherwise
      error(['MKTBOXVALID: unknown file extension ' extension]);
end; % switch extension

%%% detect ALL discontinuities, used for MKRAYDEPTH MK04122006
disconradii=mkdetectdiscon(model);


% %%% find 2740km discontinuity in IASP91 and smooth the jump out
% %%% This is to test if the peak of the time difference at 90deg is due to
% %%% this unwanted discontinuity (it is not)
% indies=find(model.z==2740);
% model.vp=model.vp([1:indies(1) (indies(2)+1):end]);
% model.vs=model.vs([1:indies(1) (indies(2)+1):end]);
% model.z=model.z([1:indies(1) (indies(2)+1):end]);

%%% Then compute the travel times
[p,a,d,deltain]=mkshootray(phase,delta,h,model);
ttcurves=mkttquick(model,h,p,a,phase,'silent');

%%% distance output of MKFINDP might contain small deviations 
%%% from desired distances. These are usually small enough to
%%% simply round them off.
rawdist=ttcurves.ttc.d;
roundeddist=round(rawdist);
ttcurves.ttc.d=roundeddist;


%%% compute deviation of computed travel times from reference
ixe=[];
ypse=[];
for indy=1:length(delta)
    indies=find(ttcurves.ttc.d==delta(indy));
    if ~isempty(indies)
       switch fullphase
          case {'PKPab'}
              minrayp=3.71; % ray parameter range fro IASP91 booklet
              maxrayp=4.44;
              miniindy=find((ttcurves.ttc.rayp(indies)>=minrayp)&(ttcurves.ttc.rayp(indies)<=maxrayp));
              comptime=ttcurves.ttc.t(indies(miniindy));
              comprayp=ttcurves.ttc.rayp(indies(miniindy));
          case {'PKPbc'}
              minrayp=2.11; % ray parameter range fro IASP91 booklet
              maxrayp=3.28;
              miniindy=find((ttcurves.ttc.rayp(indies)>=minrayp)&(ttcurves.ttc.rayp(indies)<=maxrayp));;
              comptime=ttcurves.ttc.t(indies(miniindy));
              comprayp=ttcurves.ttc.rayp(indies(miniindy));
          otherwise
              [comptime,miniindy]=min(ttcurves.ttc.t(indies));
              comprayp=ttcurves.ttc.rayp(indies(miniindy));
       end; % switch fullphase
       
       %%% catch exceptional results
       if isempty(comptime)
          comptime=NaN;
          comprayp=NaN;
       else
          %%% there should be only one value MK04122006
          comptime=comptime(1);
          comprayp=comprayp(1);
       end; % if isempty(comptime)
       
       switch takeabs
          case {'abs'}
              reftime=tt(indy);
              comptime=comptime; %ttcurves.ttc.t(indies(miniindy));
              ixe=[ixe delta(indy)];
              ypse=[ypse abs(comptime-reftime)];
              %ypse=[ypse abs(ttcurves.ttc.t(indies(miniindy))-tt(indy))];
          case {'sign'}
              reftime=tt(indy);
              comptime=comptime; %ttcurves.ttc.t(indies(miniindy));
              ixe=[ixe delta(indy)];
              ypse=[ypse comptime-reftime];
              %ypse=[ypse (ttcurves.ttc.t(indies(miniindy))-tt(indy))];
          case {'rayp'}
              reftime=tt(indy);
              comptime=comptime; %ttcurves.ttc.t(indies(miniindy));
              ixe=[ixe comprayp];
              ypse=[ypse comptime-reftime];
          case {'raydepth'}
              reftime=tt(indy);
              comptime=comptime; %ttcurves.ttc.t(indies(miniindy));
              raydep=mkraydepth(mkpdeg2rad(comprayp),...
                                model.vp,model.rp-model.z,...
                                model.rp-h,model.rp,disconradii);
              raydep=model.rp-max(raydep);
              ixe=[ixe raydep];
              ypse=[ypse comptime-reftime];
          otherwise
              error(['MKTTBOXVALID: unknown plot mode ' upper(takeabs)]);
       end; % switch takeabs
       %ypse=[ypse abs(round(100*ttcurves.ttc.t(indies(miniindy)))/100-tt(indy))];
    end; % if ~isempty(indies)
end; % for indy


%%% plot travel time deviations
figure(1);
hold on
h1=plot(ixe,ypse,linespec);
set(h1,'Tag',[fullphase ', dz=' num2str(dz) ' ' mode]);

switch takeabs
   case {'abs'}
      ylabel('|T(ttbox)-T(booklet)| [s]');
      xlabel('Epicentral Distance [deg]');
      axis auto
      ax=axis;
      axis([0 180 0 ax(4)]);
   case {'sign'}
      ylabel('T(ttbox)-T(booklet) [s]');
      xlabel('Epicentral Distance [deg]');
      axis auto
      ax=axis;
      axis([0 180 ax(3) ax(4)]);
   case {'rayp'}
      ylabel('T(ttbox)-T(booklet) [s]');
      xlabel('Ray Parameter [s/deg]');
      axis auto
   case {'raydepth'}
      ylabel('T(ttbox)-T(booklet) [s]');
      xlabel('Ray Turning Point Depth [km]');
      axis auto
   otherwise
      error(['MKTTBOXVALID: unknown plot mode ' upper(takeabs)]);
end; % switch takeabs
title(['TTBOX travel time error for model ' model.name ',  phase ' phase ', focal depth ' int2str(h) 'km.']);

box on
grid on
hold off

%%% plot travel time curves
figure(2);
hold on
h2r=plot(delta,tt,'r-');
set(h2r,'Color',[1 1 1]*0.8);
set(h2r,'LineWidth',4);
h2tt=plot(ttcurves.ttc.d,ttcurves.ttc.t,linespec);
set(h2tt,'Marker','.');
set(h2tt,'Tag',['dz=' num2str(dz) ' ' mode]);
axis auto
ax=axis;
axis([0 180 0 ax(4)]);
box on
grid on
xlabel('Epicentral Distance [deg]');
ylabel('Travel Time [s]');
title(['TTBOX travel time curves for model ' model.name  ', phase ' phase ', focal depth ' int2str(h) 'km.']);
hold off

%%% plot MKFINDP distance error
figure(3);
hold on
h3=plot(deltain,abs(deltain-d),linespec);
set(h3,'Tag',['dz=' num2str(dz) ' ' mode]);
axis auto
ax=axis;
axis([ax(1) 180 0 ax(4)]);
set(gca,'Yscale','log');
box on
grid on
xlabel('Desired Epicentral Distance [deg]');
ylabel('|Desired Minus Resulting Distance| [deg]');
title(['TTBOX shooting error for model ' model.name  ', phase ' phase ', focal depth ' int2str(h) 'km.']);
hold off



toctime=toc; disp(['MKTTBOXVALID: elapsed time: ' num2str(toctime) 'sec']);