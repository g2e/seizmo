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
h=p2(normalize(zdata1),'fgcolor','k','bgcolor','w','title',[phasename ' Aligned By PREM Arrival Time'],...
'fontsize',10,'fontweight','bold','xlabel',['Time Relative to PREM ' phasename],'ylabel','Normalized Amplitude');
saveas(h,'prealign.eps','epsc')
close(h)

% snr cut
snr=quicksnr(zdata1,nwin,swin);
goodsnr=snr(snr>minsnr);
nrecs=numel(goodsnr);
if(nrecs>2)
zdata2=taper(cut(zdata1(snr>minsnr),swin(1),swin(2)),taperwidth);
zdata3=taper(cut(envelope(zdata0(snr>minsnr)),swin(1),swin(2)),taperwidth);

% correlate
zdata4=seizmofun(zdata2,@(x)x.^recpow);
peaks=correlate(zdata4,zdata4,'npeaks',1,'absxc',false);
zdata5=seizmofun(zdata3,@(x)x.^recpow);
epeaks=correlate(zdata5,zdata5,'npeaks',1);

% correlation stats
peaks.z=fisher(peaks.cg);
epeaks.z=fisher(epeaks.cg);
peaks.zmean=nanmean(peaks.z+diag(nan(nrecs,1)));
peaks.z2std=2*sqrt(nanvariance(peaks.z+diag(nan(nrecs,1))));
peaks.rmean=ifisher(peaks.zmean);
peaks.r2stdlow=ifisher(peaks.zmean-peaks.z2std);
peaks.r2stdhigh=ifisher(peaks.zmean+peaks.z2std);
epeaks.zmean=nanmean(epeaks.z+diag(nan(nrecs,1)));
epeaks.z2std=2*sqrt(nanvariance(epeaks.z+diag(nan(nrecs,1))));
epeaks.rmean=ifisher(epeaks.zmean);
epeaks.r2stdlow=ifisher(epeaks.zmean-epeaks.z2std);
epeaks.r2stdhigh=ifisher(epeaks.zmean+epeaks.z2std);

% weights
peaks.w=peaks.z.^2.*goodsnr(:,ones(nrecs,1));
epeaks.w=epeaks.z.^2.*goodsnr(:,ones(nrecs,1));

% invert for arrivals and error (xc consistency based)
[m,Gg]=dtwalign(peaks.w,peaks.lg,1,1);
[em,eGg]=dtwalign(epeaks.w,epeaks.lg,1,1);
relarr=-m-gh(zdata2,'o');
erelarr=-em-gh(zdata2,'o');
s=dtwresid(m,peaks.w,peaks.lg,1).';
es=dtwresid(em,epeaks.w,epeaks.lg,1).';

% plot aligned
h=p2(normalize(timeshift(timeshift(zdata2,-getheader(zdata2,'o'),'io'),-relarr)),...
'fgcolor','k','bgcolor','w','title',['Aligned ' phasename],...
'fontsize',10,'fontweight','bold','ylabel','Normalized Amplitude');
saveas(h,'phase_align.eps','epsc')
close(h)
h=p2(normalize(timeshift(timeshift(zdata3,-getheader(zdata2,'o'),'io'),-erelarr)),...
'fgcolor','k','bgcolor','w','title',['Aligned ' phasename ' Envelope'],...
'fontsize',10,'fontweight','bold','ylabel','Normalized Amplitude');
saveas(h,'group_align.eps','epsc')
close(h)

% plot residuals vs dd/az
gcarc=gh(zdata2,'gcarc');
az=gh(zdata2,'az');
h=figure('color','w');
errorbar(gcarc,-m,2.*s,'--bs','MarkerEdgeColor','k',...
				'LineWidth',2,...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
hold on
errorbar(gcarc,-em,2.*es,'-g^','MarkerEdgeColor','k',...
				'LineWidth',2,...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
legend('phase','group','location','best')
title('Degree Distance vs Relative Arrival Time Residual','fontsize',10,'fontweight','bold')
ylabel('Relative Arrival Time Residual (sec)','fontsize',10,'fontweight','bold')
xlabel('Distance (degrees)','fontsize',10,'fontweight','bold')
hold off
saveas(h,'gcarc_vs_arr_time.eps','epsc')
close(h)
h=figure('color','w');
errorbar(az,-m,2.*s,'--bs','MarkerEdgeColor','k',...
				'LineWidth',2,...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
hold on
errorbar(az,-em,2.*es,'-g^','MarkerEdgeColor','k',...
				'LineWidth',2,...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
legend('phase','group','location','best')
title('Azimuth vs Relative Arrival Time Residual','fontsize',10,'fontweight','bold')
ylabel('Relative Arrival Time Residual (sec)','fontsize',10,'fontweight','bold')
xlabel('Azimuth (degrees)','fontsize',10,'fontweight','bold')
hold off
saveas(h,'az_vs_arr_time.eps','epsc')
close(h)

% plot amplitudes vs dd/az
phamp=gh(zdata2,'depmax')-gh(zdata2,'depmin');
envamp=gh(zdata3,'depmax');
h=figure('color','w');
errorbar(gcarc,phamp,phamp./goodsnr,'--bs',...
				'LineWidth',2,...
				'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
hold on
errorbar(gcarc,envamp,envamp./goodsnr,'-g^',...
				'LineWidth',2,...
				'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
set(gca,'yscale','log');
legend('phase','group','location','best')
title('Degree Distance vs Amplitude','fontsize',10,'fontweight','bold')
ylabel('Amplitude (nm)','fontsize',10,'fontweight','bold')
xlabel('Distance (degrees)','fontsize',10,'fontweight','bold')
hold off
saveas(h,'gcarc_vs_log_amp.eps','epsc')
close(h)
h=figure('color','w');
errorbar(az,phamp,phamp./goodsnr,'--bs',...
				'LineWidth',2,...
				'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
hold on
errorbar(az,envamp,envamp./goodsnr,'-g^',...
				'LineWidth',2,...
				'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
set(gca,'yscale','log');
legend('phase','group','location','best')
title('Azimuth vs Amplitude','fontsize',10,'fontweight','bold')
ylabel('Amplitude (nm)','fontsize',10,'fontweight','bold')
xlabel('Azimuth (degrees)','fontsize',10,'fontweight','bold')
hold off
saveas(h,'az_vs_log_amp.eps','epsc')
close(h)

% output
rep=ones(numel(zdata2),1);
stidx=gh(zdata2,'resp0');
blah=[evidx(rep,1) stidx bankidx(rep,1) bandidx(rep,1) arridx(rep,1) nwin(rep,:) swin(rep,:) taperwidth(rep,1) recpow(rep,1) goodsnr peaks.rmean(:) m relarr s phamp epeaks.rmean(:) em erelarr es envamp];
save([event '_' phasename '_f' num2str(round(sp)) '-' num2str(round(lp)) '_w' num2str(swin(1)) '-' num2str(swin(2)) '_snr' num2str(minsnr) '.align'],'blah','-ascii');
blah
end

