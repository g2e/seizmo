function [psdgram]=noise_psdgram(indir,length,overlap,func)
%NOISE_PSDGRAM    Make power spectra for noise analysis
%
%    Usage:    psdgram=noise_psdgram(indir)
%              psdgram=noise_psdgram(indir,length,overlap,func)
%
%    Description:
%     PSDGRAM=NOISE_PSDGRAM(INDIR) creates power spectral density (PSD)
%     spectrograms for each component in the noise directory INDIR.  INDIR
%     should be the output from NOISE_SETUP (i.e., the records should be in
%     displacement).  The output PSDGRAM is a struct with one element per
%     component in INDIR and uses the format of PLOTPSDGRAM which is
%     typically used to analyze ocean wave spectra.  The output is also
%     written within INDIR as INDIR/noise_psdgrams_DATETIMESTRING.mat.
%
%     PSDGRAM=NOISE_PSDGRAM(INDIR,LENGTH,OVERLAP,FUNC) allows specifying
%     the length and overlap of the PSDs used to compute the spectra for
%     each timesection.  This basically allows you to use sub-windows for
%     the spectral analysis.  LENGTH and OVERLAP are in minutes.  FUNC is
%     the function used to compute the timesection spectra from the
%     sub-window spectra and should operate on each column of the input
%     matrix separately (1 column per frequency).  The default LENGTH is
%     the entire timesection length (can be specified with []).  The
%     default OVERLAP is 0.  The default FUNC is @median (@mean, @min, &
%     @max are possible alternatives).
%
%    Notes:
%
%    Examples:
%     % The default from NOISE_SETUP is 180 minute timesections with no
%     % overlap.  Setting a psd length of 10 minutes (no overlap) gives
%     % 18 sub-windows (out of which the median value for each frequency
%     % is used to make the output spectra).  The spectra are then written
%     % to the input directory:
%     noise_psdgram('my_setup_dir',10);
%
%    See also: PLOTPSDGRAM, CHKPSDGRAM, READ_NDBC_SWDEN, NOISE_SETUP

%     Version History:
%        Apr. 15, 2013 - initial version, no 50% data in subwindow bugfix
%        Apr. 17, 2013 - use powerspectraldensity to get psd
%        Jan. 26, 2014 - abs path exist fix
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% directory separator
fs=filesep;

% check input directory
if(~isstring(indir))
    error('seizmo:noise_psdgram:fileNotString',...
        'INDIR must be a string!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~exist(indir,'dir'))
    error('seizmo:noise_psdgram:dirConflict',...
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

% convert time ranges to numeric arrays
t=char(tsdirs);
tsbgn=str2double([cellstr(t(:,1:4)) cellstr(t(:,6:8)) ...
    cellstr(t(:,10:11)) cellstr(t(:,13:14)) cellstr(t(:,16:17))]);
tsend=str2double([cellstr(t(:,19:22)) cellstr(t(:,24:26)) ...
    cellstr(t(:,28:29)) cellstr(t(:,31:32)) cellstr(t(:,34:35))]);
tsmid=(gregorian2serial(tsbgn)+gregorian2serial(tsend))/2;
tslen=timediff(tsbgn(1,:),tsend(1,:));
clear t;

% make sure all have the same amount of time
if(numel(unique(timediff(tsbgn,tsend)))~=1)
    error('seizmo:noise_psdgram:corruptDir',...
        'Time section directories under INDIR are not consistent!');
end

% default length, overlap, func
if(nargin<2 || isempty(length)); length=tslen; end
if(nargin<3 || isempty(overlap)); overlap=0; end
if(nargin<4 || isempty(func)); func=@median; end

% check length, overlap, func
if(~isnumeric(length) || ~isreal(length) ...
        || ~isscalar(length) || length~=fix(length))
    error('seizmo:noise_psdgram:badInput',...
        'LENGTH must be a real-valued scalar integer in minutes!');
elseif(length<=0 || length>tslen)
    error('seizmo:noise_psdgram:badInput',...
        'LENGTH must be >0 and no more than the timesection length!');
elseif(~isnumeric(overlap) || ~isreal(overlap) ...
        || ~isscalar(overlap) || overlap~=fix(overlap))
    error('seizmo:noise_psdgram:badInput',...
        'OVERLAP must be a real-valued scalar integer in minutes!');
elseif(overlap<0 || overlap>=tslen)
    error('seizmo:noise_psdgram:badInput',...
        'OVERLAP must be less than the timesection length and >=0!');
elseif(~isa(func,'function_handle') || ~isscalar(func))
    error('seizmo:noise_psdgram:badInput',...
        'FUNC must be a function handle!');
end

% detail message
n=numel(tsdirs);
verbose=seizmoverbose(false);
if(verbose)
    disp('MAKING PSD SPECTROGRAMS');
    print_time_left(0,n);
end

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% loop over tsdirs
sname=cell(0,1);
for i=1:n
    % read in timesection data
    try
        try
            data=load([indir fs tsdirs{i}(1:4) fs tsdirs{i} ...
                fs 'noise_records.mat']);
            data=data.noise_records;
        catch
            data=readseizmo([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
        end
    catch
        continue;
    end
    if(isempty(data)); continue; end
    
    % check the samplerate
    [dt,knetwk,kstnm,khole,kcmpnm]=getheader(data,...
        'delta','knetwk','kstnm','khole','kcmpnm');
    if(numel(unique(dt))~=1)
        error('seizmo:noise_psdgram:badData',...
            'Not all records have the same samplerate!');
    end
    if(~exist('dt0','var'))
        dt0=dt(1);
        swlenmax=ceil(tslen/dt0);
        swlen=ceil(length*60/dt0);
        swstep=round(swlen-ceil(overlap*60/dt0));
        nsw=floor((swlenmax-swlen)/swstep+1);
        swlen2=2^nextpow2(swlen);
    elseif(dt0~=dt(1))
        error('seizmo:noise_psdgram:badData',...
            'Sample rate of records changes between timesections!');
    end
    
    % get knames
    kname=lower(strcat(knetwk,'.',kstnm,'.',khole,'.',kcmpnm));
    
    % get the spectra
    swspectra=cell(numel(data),1);
    for j=1:nsw
        % cut
        b=(j-1)*60*(length-overlap);
        [swdata,failed]=cut(data,b,b+length*60); % assumes tsbgn==0
        ok=true(numel(data),1);
        ok(failed)=false;
        oki=find(ok);
        
        % handle none
        if(isempty(swdata)); continue; end
        
        % require 50% of subwindow
        npts=getheader(swdata,'npts');
        failed=npts<swlen/2;
        oki(failed)=[];
        swdata(failed)=[];
        
        % handle none
        if(isempty(swdata)); continue; end
        
        % remove trend & taper
        swdata=taper(removetrend(swdata),.5);
        
        % pad with zeros
        swdata=cut(swdata,'x',1,'n',swlen2,'fill',true);
        
        % get power spectra (not in dBs)
        swdata=solofun(powerspectraldensity(divide(swdata,1e9)),...
            @(x)10.^(x/10));
        
        % grab frequency once
        if(~exist('f','var'))
            [npts,delta]=getheader(swdata(1),'npts','delta');
            f=(0:delta:delta*(npts-1)).';
        end
        
        % append spectra to storage cells (one per record)
        for k=1:numel(oki)
            swspectra{oki(k)}=[swspectra{oki(k)} swdata(k).dep];
        end
    end
    
    % apply function to each cell
    % - avoid call for single spectra to avoid issues & for speed
    % - transpose so freq down columns
    ok=true(numel(data),1);
    for j=1:numel(data)
        if(isempty(swspectra{j})); ok(j)=false; continue; end
        swspectra{j}=swspectra{j}';
        if(size(swspectra{j},2)==1); continue; end
        swspectra{j}=func(swspectra{j});
    end
    
    % remove unused knames & spectra
    kname(~ok)=[];
    swspectra(~ok)=[];
    
    % output result to another cell array based on names
    % depending on knames
    % 1. append spectra
    % 2. create new spectra
    [tf,loc]=ismember(kname,sname);
    if(any(tf))
        loc=loc(tf);
        tf=find(tf);
        for j=1:numel(tf)
            psdgram(loc(j)).spectra=[psdgram(loc(j)).spectra;
                swspectra{tf(j)}];
            psdgram(loc(j)).time=[psdgram(loc(j)).time; tsmid(i)];
        end
    end
    
    % create new spectra
    if(any(~tf))
        npsd=numel(sname);
        sname=[sname; kname(~tf)];
        new=find(~tf);
        for j=1:numel(new)
            psdgram(npsd+j).name=kname{new(j)};
            psdgram(npsd+j).units='m^2/Hz';
            psdgram(npsd+j).time=tsmid(i);
            psdgram(npsd+j).freq=f;
            psdgram(npsd+j).spectra=swspectra{new(j)};
        end
    end
    
    % detail message
    if(verbose); print_time_left(i,n); end
end

% toggle checking back
seizmocheck_state(oldseizmocheckstate);
checkheader_state(oldcheckheaderstate);

% reset verbosity
seizmoverbose(verbose);

% save
save([indir fs 'noise_psdgrams_' datestr(now,30) '.mat'],'psdgram');

end
