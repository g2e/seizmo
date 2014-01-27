function []=noise_mat2sac(noisedir,rmflag)
%NOISE_MAT2SAC    Convert directories of noise data from mat to sac files
%
%    Usage:    noise_mat2sac(noisedir,rmflag)
%
%    Description:
%     NOISE_MAT2SAC(NOISEDIR,RMFLAG) converts the noise_records.mat files
%     within the noise directory structure in NOISEDIR to sac files within
%     the same directory location.  The directories are expected to be
%     formatted as for NOISE_SETUP & NOISE_PROCESS.  The logical RMFLAG
%     indicates if the noise_records.mat files are to be removed as well -
%     the default value (true) removes them.
%
%    Notes:
%
%    Examples:
%     % The noise_records.mat files are not very portable, use
%     % noise_mat2sac to convert them to SAC files instead (here we do not
%     % delete .mat files as we are only inspecting):
%     noise_mat2sac('setup-3hr',false)
%
%    See also: NOISE_SAC2MAT, NOISE_SETUP, NOISE_PROCESS, NOISE_STACK

%     Version History:
%        Mar. 26, 2013 - initial version
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 26, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% directory separator
fs=filesep;

% check directory
if(~isstring(noisedir))
    error('seizmo:noise_mat2sac:fileNotString',...
        'NOISEDIR must be a string!');
end
if(~isabspath(noisedir)); noisedir=[pwd fs noisedir]; end
if(~exist(noisedir,'dir'))
    error('seizmo:noise_mat2sac:dirConflict',...
        ['Noise Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],noisedir);
end

% default/check remove flag
if(nargin<2 || isempty(rmflag)); rmflag=true; end
if(~isscalar(rmflag) || ~islogical(rmflag))
    error('seizmo:noise_mat2sac:dirConflict',...
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
    disp('Converting noise_records.mat to SAC files');
    print_time_left(0,ndirs);
end

% read, delete, write
for i=1:numel(tsdirs)
    try
        noise_records=load(strcat(noisedir,fs,tsdirs{i}(1:4),fs,...
            tsdirs{i},fs,'noise_records.mat'));
        noise_records=noise_records.noise_records;
    catch
        continue;
    end
    writeseizmo(noise_records,'path',...
        strcat(noisedir,fs,tsdirs{i}(1:4),fs,tsdirs{i},fs));
    if(rmflag)
        delete(strcat(noisedir,fs,tsdirs{i}(1:4),fs,...
            tsdirs{i},fs,'noise_records.mat'));
    end
    if(verbose); print_time_left(i,ndirs); end
end

% return verbosity
seizmoverbose(verbose);

end
