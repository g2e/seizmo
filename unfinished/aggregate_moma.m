function []=aggregate_moma(event,phase,n)
% for slowdecay processing

% event to process
%event=1;
%phase='Pdiff';

% load info
sta=load('moma.station');
evt=load('moma.event');
evrow=event==evt(:,1);
evla=evt(evrow,7);
evlo=evt(evrow,8);
evdir=[num2str(evt(evrow,2)-1900) num2str(evt(evrow,3),'%03d')];

% filter info
period=1./[
   1.0000000e-02           8.0000000e-03           1.2000000e-02
   1.2000000e-02           9.6000000e-03           1.4400000e-02
   1.4400000e-02           1.1520000e-02           1.7280000e-02
   1.7280000e-02           1.3824000e-02           2.0736000e-02
   2.0736000e-02           1.6588800e-02           2.4883200e-02
   2.4883200e-02           1.9906560e-02           2.9859840e-02
   2.9859840e-02           2.3887872e-02           3.5831808e-02
   3.5831808e-02           2.8665446e-02           4.2998170e-02
   4.2998170e-02           3.4398536e-02           5.1597804e-02
   5.1597804e-02           4.1278243e-02           6.1917364e-02
   6.1917364e-02           4.9533891e-02           7.4300837e-02
   7.4300837e-02           5.9440670e-02           8.9161004e-02
   8.9161004e-02           7.1328804e-02           1.0699321e-01
   1.0699321e-01           8.5594564e-02           1.2839185e-01
   1.2839185e-01           1.0271348e-01           1.5407022e-01
   1.5407022e-01           1.2325617e-01           1.8488426e-01
   1.8488426e-01           1.4790741e-01           2.2186111e-01];

% get station distances
gcarc=sphericalinv(evla,evlo,sta(:,2),sta(:,3));

% enter event directory
cd(evdir)

% enter processing directory
cd([phase '.normal'])

% load results
a=cell(17,1); fbad=false(1,17);
for i=1:17
    si=sprintf('%02d',i);
    file=xdir([si '*' filesep '*.align']);
    fbad(i)=numel(file)~=1;
    if(~fbad(i))
        file=[file.path filesep file.name];
        a{i}=load(file);
    end
end

% exit processing directory
cd ..

% make dd vs tt
for i=find(~fbad)
    idx{i}=a{i}(:,2)+1;
    tt{i}=a{i}(:,15);
    amp{i}=a{i}(:,17);
    dist{i}=gcarc(idx{i});
end

% read in corrections
switch event
    case 1
        cidx=[15 6 2 17 9 12 18 5 10 11 7 8 19 0]+1;
    case 2
        cidx=[17 5 6 10 13 2 4 14 11 9 1 8 12 16 18 7 3 19 0]+1;
    case 3
        cidx=[5 3 6 11 10 2 9 16 8 7 17 1 4 18 19 0]+1;
    case 4
        cidx=[7 3 1 2 9 8 14 16 17 4 18 5 10 0]+1;
    case 5
        cidx=[3 2 8 9 12 14 11 6 17 7 10 13 4 18 1 16 5 0]+1;
    case 6
        cidx=[7 11 5 3 17 4 6 2 9 15 18 13 12 16 8 14 0]+1;
    case 7
        cidx=[0 6 5 13 16 9 12 18 4 15 11 7 3 17 8 14]+1;
end
switch phase
    case 'Pdiff'
        mantlefile='Pdiff.mantle.hmsl-p06.upswing';
        crustfile='Pdiff.crust';
        elevfile='Pdiff.elev';
        ellipfile='Pdiff.ellip';
    case {'SHdiff' 'SVdiff'}
        mantlefile='Sdiff.mantle.hmsl-s06.upswing';
        crustfile='Sdiff.prem.crust';
        elevfile='Sdiff.prem.elev';
        ellipfile='Sdiff.ellip.ak135';
end
txt=readtxt(mantlefile);
words=reshape(getwords(txt),2,[])';
corr(cidx,1)=str2double(words(:,2));
txt=readtxt(crustfile);
words=reshape(getwords(txt),2,[])';
corr(cidx,2)=str2double(words(:,2));
txt=readtxt(elevfile);
words=reshape(getwords(txt),2,[])';
corr(cidx,3)=str2double(words(:,2));
txt=readtxt(ellipfile);
words=reshape(getwords(txt),2,[])';
corr(cidx,4)=str2double(words(:,2));

% get spreading correction
gs=sin(gcarc*pi/180).^-0.5;

% add corrections
ctt=cell(17,1); camp=ctt;
slow=nan(17,1); decay=nan(17,1);
for i=find(~fbad)
    % add corrections
    ctt{i}=tt{i}-sum(corr(idx{i},:),2);
    camp{i}=amp{i}./gs(idx{i});
    
    % get slowness
    [dist{i} amp{i} gs(idx{i}) amp{i}./gs(idx{i})];
    [dist{i} tt{i} sum(corr(idx{i},:),2) ctt{i}];
    tmp=wlinem(dist{i},ctt{i});
    slow(i)=tmp(2);
    
    % get closest station
    [tmp,minidx]=min(dist{i});
    mamp=camp{i}(minidx);
    
    % get decays
    log(camp{i}/mamp);
    tmp=wlinem(dist{i},log(camp{i}/mamp));
    decay(i)=-tmp(2);
end
slow
decay

% make plots
figure;
plot(period(n,1),slow(n),'g','marker','^','markerfacecolor','b',...
    'markersize',10,'markeredgecolor','g','linewidth',2);
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold');
grid on
xlabel('Period (s)','fontsize',12,'fontweight','bold');
ylabel('Slowness (s/deg)','fontsize',12,'fontweight','bold');
title([evdir ' - ' phase ' - Slowness'],...
    'fontsize',12,'fontweight','bold');
switch phase
    case 'Pdiff'
        ylim([4.2 5.2]);
    case {'SHdiff' 'SVdiff'}
        ylim([8.2 9.2]);
end
xlim([0 100]);

figure;
plot(period(n,1),decay(n),'g','marker','s','markerfacecolor','b',...
    'markersize',10,'markeredgecolor','g','linewidth',2);
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold');
grid on
xlabel('Period (s)','fontsize',12,'fontweight','bold');
ylabel('Decay Constant','fontsize',12,'fontweight','bold');
title([evdir ' - ' phase ' - Amplitude Decay'],...
    'fontsize',12,'fontweight','bold');
ylim([0 0.15]);
xlim([0 100]);

% exit event directory
cd ..

end

