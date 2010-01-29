% parameters (set by calling script/function)
%evidx=2;		% event idx (chrono)
%bankidx=3;		% filter bank index (vw, w, n, s)
%bandidx=3;		% filter band index
%arridx=1;		% arrival index (Pdiff=1, SHdiff=2, SVdiff=3)
%owin=[-200 400];% window to decide following windows
%nwin=[-200 50]; % noise window
%swin=[50 300];	% signal window
%taperwidth=0.2; % taper width
%minsnr=3; 		% snr must be > minsnr
%recpow=3; 		% multiply recs by this before correlating, must be positive and odd

% get event name
events={'95078' '95136' '95175' '95180' '95184' '95235' '96048'};
event=events{evidx};

% get phase name / component name
phases={'Pdiff' 'Sdiff' 'Sdiff'};
phasenames={'Pdiff' 'SHdiff' 'SVdiff'};
phase=phases{arridx};
phasename=phasenames{arridx};
cmps={'BHZ' 'BHT' 'BHR'};
cmp=cmps{arridx};

% get filter
banklimits=[0.01 0.2];
bankoption='variable';
bankparams=[1.6 0.8; 0.8 0.4; 0.4 0.2; 0.2 0.1];
bank=1./filter_bank(banklimits,bankoption,bankparams(bankidx,1),bankparams(bankidx,2));
cp=bank(bandidx,1);
lp=bank(bandidx,2);
sp=bank(bandidx,3);

% initial window
cd(['../' sprintf('%02d',bandidx) '_CENTER_' num2str(cp) '_LOWPASS_' num2str(sp) '_HIGHPASS_' num2str(lp)]);
zdata=r(['*' cmp '*']);
zdata=sortbyfield(zdata,'gcarc');
[arrtime,idx]=getarrival(zdata,phase);
zdata0=timeshift(zdata,-arrtime,strcat({'it'},num2str(idx)));
zdata1=cut(zdata0,owin(1),owin(2));
[zdata2,win,fh]=userwindow(normalize(zdata1),true,@removemean,'fgcolor','k','bgcolor','w','title',[phasename ' Aligned By PREM Arrival Time'],...
'fontsize',10,'fontweight','bold','xlabel',['Time Relative to PREM ' phasename],'ylabel','Normalized Amplitude');
saveas(fh(1),'prealign.eps','epsc')
saveas(fh(2),'windowed.eps','epsc')
win=win(1,:)
save window.mat win -ascii
close(fh)

