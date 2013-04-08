function [psdgram]=noise_psdgram(indir)

% get year directories and time-section directories
fs=filesep;
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
yrdir={dirs.name};
nyr=numel(yrdir);
tsdir=cell(size(yrdir));
for i=1:nyr
    dirs=xdir([indir fs yrdir{i}]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    tsdir{i}={dirs.name};
end
tsdirs=[tsdir{:}]'; % all timesections
clear dirs yrdir nyr tsdir;

% convert time ranges to numeric arrays
t=char(tsdirs);
tsbgn=str2double([cellstr(t(:,1:4)) cellstr(t(:,6:8)) ...
    cellstr(t(:,10:11)) cellstr(t(:,13:14)) cellstr(t(:,16:17))]);
tsend=str2double([cellstr(t(:,19:22)) cellstr(t(:,24:26)) ...
    cellstr(t(:,28:29)) cellstr(t(:,31:32)) cellstr(t(:,34:35))]);
clear t;

psdgram.time=(gregorian2serial(tsbgn)+gregorian2serial(tsend))/2;
psdgram.freq=[];
psdgram.spectra=[];
psdgram.units='m^2/Hz';
n=numel(tsdirs);
print_time_left(0,n);
verbose=seizmoverbose(false);
% loop over tsdirs, fft, stack psd
for i=1:n
    data=load(strcat(indir,fs,tsdirs{i}(1:4),fs,...
        tsdirs{i},fs,'noise_records.mat'));
    data=data.noise_records;
    data=divide(data,1e9);
    data=divide(addrecords(keepam(dft(taper(removetrend(...
        cut(data,'x',1,'n',8192,'fill',true)),.5)))),numel(data));
    if(i==1)
        [b,npts,delta]=getheader(data,'b','npts','delta');
        psdgram.freq=(0:npts-1)'*delta;
    end
    psdgram.spectra(i,:)=data.dep;
    print_time_left(i,n);
end
seizmoverbose(verbose);

% save
save(fullfile(indir,'noise_psd.mat'),'psdgram');

end
