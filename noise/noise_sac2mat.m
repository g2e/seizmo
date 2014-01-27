function []=noise_sac2mat(noisedir,rmflag)
%NOISE_SAC2MAT    Convert directories of noise data from sac to mat files
%
%    Usage:    noise_sac2mat(noisedir,rmflag)
%
%    Description:
%     NOISE_SAC2MAT(NOISEDIR,RMFLAG) converts the sac files within the
%     noise directory structure in NOISEDIR to noise_records.mat files
%     within the same directory location.  The directories are expected to
%     be formatted as for NOISE_SETUP & NOISE_PROCESS.  The logical RMFLAG
%     indicates if the sac files are to be removed as well - the default
%     value (true) removes them.
%
%    Notes:
%     - This is a fairly dumb script that will just read in everything it
%       can from each time window directory, remove and write.  Use at your
%       own risk!
%
%    Examples:
%     % Sometimes there are too many sac files and it starts to hurt the
%     % filesystem itself due to fragmentation.  Converting to mat files
%     % packages all the sac files for each time window into 1 file:
%     noise_sac2mat('setup-3hr')
%
%    See also: NOISE_MAT2SAC, NOISE_SETUP, NOISE_PROCESS, NOISE_STACK

%     Version History:
%        Mar. 26, 2013 - initial version
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% directory separator
fs=filesep;

% check directory
if(~isstring(noisedir))
    error('seizmo:noise_sac2mat:fileNotString',...
        'NOISEDIR must be a string!');
end
if(~isabspath(noisedir)); noisedir=[pwd fs noisedir]; end
if(~exist(noisedir,'dir'))
    error('seizmo:noise_sac2mat:dirConflict',...
        ['Noise Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end

% default/check remove flag
if(nargin<2 || isempty(rmflag)); rmflag=true; end
if(~isscalar(rmflag) || ~islogical(rmflag))
    error('seizmo:noise_sac2mat:dirConflict',...
        'RMFLAG must be TRUE or FALSE!');
end

% get year directories and time-section directories
dirs=xdir([noisedir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
yrdir={dirs.name};
nyr=numel(yrdir);
tsdir=cell(size(yrdir));
for i=1:nyr
    dirs=xdir([noisedir fs yrdir{i}]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    tsdir{i}={dirs.name};
end
tsdirs=[tsdir{:}]'; % all timesections
clear dirs yrdir nyr tsdir;

% progress bar
verbose=seizmoverbose(false);
if(verbose)
    ndirs=numel(tsdirs);
    disp('Converting SAC files to noise_records.mat');
    print_time_left(0,ndirs);
end

% read, delete, write
for i=1:numel(tsdirs)
    try
        noise_records=readseizmo(...
            strcat(noisedir,fs,tsdirs{i}(1:4),fs,tsdirs{i}));
    catch
        continue;
    end
    save(strcat(noisedir,fs,tsdirs{i}(1:4),fs,tsdirs{i},fs,...
        'noise_records.mat'),'noise_records');
    if(rmflag)
        for j=1:numel(noise_records)
            delete(strcat(noisedir,fs,tsdirs{i}(1:4),fs,...
                tsdirs{i},fs,noise_records(j).name));
        end
    end
    if(verbose); print_time_left(i,ndirs); end
end

% return verbosity
seizmoverbose(verbose);

end
