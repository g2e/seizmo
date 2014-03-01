function []=noise_rms_calc(indir)
%NOISE_RMS_CALC    Calculates RMS values for noise analysis
%
%    Usage:    noise_rms_calc(indir)
%
%    Description:
%     NOISE_RMS_CALC(INDIR) computes the root-mean-square (rms) for each
%     record in every timesection under noise analysis directory structure
%     INDIR.  See NOISE_SETUP for more info on the layout of INDIR.  The
%     rms info is saved as noise_rms_info.mat under the directory INDIR.
%
%    Notes:
%
%    Examples:
%     % Lookout for bad data early in the process:
%     noise_setup('raw','3hr')
%     noise_rms_calc('3hr')
%     noise_process('3hr','ncfs',[],'rm','rms','rc',5);
%
%    See also: NOISE_SETUP, NOISE_PROCESS, NOISE_STACK, NOISE_OVERVIEW

%     Version History:
%        Aug. 27, 2012 - initial version
%        Aug. 30, 2012 - output is now a seizmo dataset
%        Aug. 31, 2012 - fixed to handle time gaps
%        July 24, 2013 - update .path field to '.'
%        Aug.  8, 2013 - fixed buggy time step code using function GCDN,
%                        tsbgn & tsend output now include times for every
%                        point in the output data for sanity
%        Sep. 24, 2013 - fixed example to avoid xc directory issue
%        Jan. 26, 2014 - abs path exist fix
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% directory separator
fs=filesep;

% check directory
if(~isstring(indir))
    error('seizmo:noise_rms_calc:fileNotString',...
        'INDIR must be a string!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~exist(indir,'dir'))
    error('seizmo:noise_rms_calc:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end

% get year directories and time-section directories
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
nts=numel(tsdirs);
clear dirs yrdir nyr tsdir;

% convert time ranges to numeric arrays
t=char(tsdirs);
tsbgn=gregorian2serial(str2double([cellstr(t(:,1:4)) cellstr(t(:,6:8)) ...
    cellstr(t(:,10:11)) cellstr(t(:,13:14)) cellstr(t(:,16:17))]));
tsend=gregorian2serial(str2double([cellstr(t(:,19:22)) cellstr(t(:,24:26)) ...
    cellstr(t(:,28:29)) cellstr(t(:,31:32)) cellstr(t(:,34:35))]));
clear t;

% check timesection lengths
% - adding UTC to timediff here would break the detection below
if(numel(unique(tsend-tsbgn))>1)
    error('seizmo:noise_rms_calc:variableTimeSectionWidth',...
        'Timesection time spans vary in INDIR!');
end
tslen=round((tsend(1)-tsbgn(1))*86400);

% need to get time step
tsstep=gcdn(round(1440*diff(tsbgn)))*60;
tsidx=round((tsbgn-tsbgn(1))*86400)/tsstep+1;

% temporary verbose control
verbose=seizmoverbose(false);
if(verbose)
    disp('Calculating RMS Values');
    print_time_left(0,nts);
end

% loop over time section directories
name=cell(0,1);
rms=nan(nts,0);
ndata=0;
for i=1:nts
    % read in timesection data
    try
        tsdata=readseizmo(strcat(indir,fs,tsdirs{i}(1:4),fs,tsdirs{i},fs));
    catch
        % no data...
        rms(i,:)=nan;
        continue;
    end
    
    % get names
    [nidx,kname]=getcomponentidx(tsdata);
    
    % debug
    if(numel(kname)~=numel(tsdata))
        disp('Multiple files for same KNAME in timesection!');
        disp(tsdirs{i});
        redraw=true;
    else
        redraw=false;
    end
    
    % find name indexing
    [tf,idx]=ismember(kname,name);
    
    % address issues for new stations
    if(any(~tf))
        idx(~tf)=ndata+(1:sum(~tf));
        rms(1:i-1,idx(~tf))=nan;
        name=[name; kname(~tf)]; %#ok<AGROW>
    end
    
    % calc rms
    rms(i,:)=nan;
    rms(i,idx(nidx))=getvaluefun(tsdata,@(x)sqrt(mean((x-mean(x)).^2)));
    
    % data structure for output
    if(any(~tf))
        [tsdata.dep]=deal([]);
        if(ndata==0)
            data(idx(nidx))=tsdata;
        else
            % the indexing is insane here...
            data(idx(nidx(~tf(nidx))))=tsdata(~tf(nidx));
        end
        ndata=numel(data);
    end
    
    % detail message
    if(verbose); print_time_left(i,nts,redraw); end
end

% insert rms & timing info into data
npts=max(tsidx);
[depmin,depmen,depmax]=deal(nan(ndata,1));
for i=1:ndata
    data(i).dep=nan(npts,1);
    data(i).dep(tsidx)=rms(:,i);
    depmin(i)=min(data(i).dep);
    depmen(i)=nanmean(data(i).dep);
    depmax(i)=max(data(i).dep);
end
data=changeheader(data,'z6',serial2gregorian(tsbgn(1)),...
    'depmin',depmin,'depmen',depmen,'depmax',depmax,...
    'b',0,'delta',tsstep,'npts',npts,'e',(npts-1)*tsstep,'f',tslen);
data=changename(data,'name',name);
data=changepath(data,'path','.');

% update tsbgn & tsend to include skipped (no data) time windows
tsbgn=tsbgn(1)+(0:tsstep/86400:tsstep/86400*(npts-1));
tsend=tsbgn+tslen/86400;

% return verbosity
seizmoverbose(verbose);

% save output
save([indir fs 'noise_rms_info.mat'],'data','tsbgn','tsend');

end
