function [ok]=webinstall_mmap(mypath)
%WEBINSTALL_MMAP    Install M_Map components
%
%    Usage:    ok=webinstall_mmap
%              ok=webinstall_mmap(path)
%
%    Description:
%     OK=WEBINSTALL_MMAP downloads the M_Map zip file to the directory
%     where this mfile is located, extracts its contents to the 'm_map'
%     directory and adds/saves it to the path.  The download is only 0.6
%     megabytes.  Just be aware that you will not be able to do anything
%     else at the Matlab/Octave prompt while waiting for the file to
%     download.
%
%     OK=WEBINSTALL_MMAP(PATH) installs the M_Map toolbox under the
%     directory given by PATH.
%
%    Notes:
%
%    Examples:
%     % Update M_Map by uninstalling, deleting the zip file & reinstalling:
%     uninstall_mmap;
%     delete([fileparts(which('webinstall_mmap')) filesep 'm_map1.4.zip']);
%     webinstall_mmap
%
%    See also: UNINSTALL_MMAP, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXPORTFIG, WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS,
%              WEBINSTALL_IRISWS, UNINSTALL_IRISWS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO, MMAP,
%              M_MAP

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - add m_map_fixes, doc update, flip savepath logic
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Apr. 25, 2012 - copy m_map/private to m_map_fixes/private
%        Jan. 14, 2014 - throw warning if problem
%        Jan. 15, 2014 - updated example to actually update M_Map, updated
%                        See also list
%        Jan. 27, 2014 - added isabspath for abs path fix to path option
%        Feb. 27, 2014 - updated See also list
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% directory separator
fs=filesep;

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=fileparts(mfilename('fullpath'));
else
    if(~isabspath(mypath)); mypath=[pwd fs mypath]; end
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_mmap:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt M_Map install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current version
    mmap='m_map1.4.zip';
    
    % grab file
    url='http://www.eos.ubc.ca/%7Erich/';
    disp([' Getting ' mmap]);
    if(exist(mmap,'file'))
        if(~exist([mypath fs mmap],'file'))
            copyfile(which(mmap),'.');
        end
    else
        urlwrite([url mmap],mmap);
    end
    
    % delete pre-existing directory if in Octave
    mmapdir=[mypath fs 'm_map'];
    if(exist(mmapdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',mmapdir);
        y=rmdir(mmapdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install M_Map
    unzip(mmap);
    addpath(mmapdir);
    if(exist([mypath fs 'm_map_fixes'],'dir'))
        addpath([mypath fs 'm_map_fixes']);
        copyfile([mmapdir fs 'private'],...
            [mypath fs 'm_map_fixes' fs 'private']);
    end
    savepath;
    
    % return
    cd(cwd);
    
    % all good
    ok=true;
catch
    le=lasterror;
    warning(le.identifier,le.message);
    ok=false;
    cd(cwd);
end

end

function [lgc]=isabspath(path,iswindows)
%ISABSPATH    Determines if a path is an absolute path or not
%
%    Usage:    lgc=isabspath(path)
%              lgc=isabspath(path,iswindows)
%
%    Description:
%     LGC=ISABSPATH(PATH) checks if the path(s) in PATH are relative or
%     absolute and returns TRUE for those that are absolute paths.  PATHS
%     may be a string, char array or a cell string array.  LGC is a logical
%     array with one element per path in PATH.  This is useful to find
%     relative paths so you can convert them to absolute paths for
%     functions like EXIST.  The determination is done by discovering the
%     OS type of the current system using ISPC.
%
%     LGC=ISABSPATH(PATH,ISWINDOWS) allows setting the OS type for
%     determining if the paths are absolute or not when the paths are not
%     valid paths for the current machine.  For instance, set ISWINDOWS to
%     FALSE for Unix, Linux or MACOSX paths when you are using MicrosoftTM
%     WindowsTM.  ISWINDOWS must be TRUE or FALSE (scalar only).
%
%    Notes:
%     - The path is not required to exist or even to be valid!  This just
%       does a simple test on each path given the OS (e.g., is the first
%       character a '/' for unix).
%
%    Examples:
%     % Test a few relative paths:
%     isabspath('./somedir')
%     isabspath('../somedir')
%     isabspath('~/somedir')
%     isabspath('..\somewindir')
%
%     % Test a few absolute paths:
%     isabspath('/home')
%     isabspath('/usr/share/../bin')
%     isabspath('c:\Programs')
%
%     % And a few invalid ones:
%     isabspath('/\')                   % absolute path to the '\' dir?
%     isabspath('somedir\c:/somewhere') % win drive in a unix dir under pwd
%     isabspath('\\someserver\somedir') % maybe you can add this feature...
%
%    See also: ISPC, ISUNIX

%     Version History:
%        Jan. 27, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 11:15 GMT

% todo:

% check number of inputs
error(nargchk(1,2,nargin));

% check/fix path
if(ischar(path))
    path=cellstr(path);
elseif(~iscellstr(path))
    error('seizmo:isabspath:badInput',...
        'PATH must be a string, char array or a cellstr array!');
end

% check/default os
if(nargin<2 || isempty(iswindows)); iswindows=ispc; end
if(~islogical(iswindows) || ~isscalar(iswindows))
    error('seizmo:isabspath:badInput',...
        'ISWINDOWS must be TRUE or FALSE!');
end

% preallocate output as all relative paths
lgc=false(size(path));

% act by os
if(iswindows) % windows
    for i=1:numel(path)
        % require drive char to be a-z,A-Z
        if(isempty(path{i})); continue; end
        drive=double(upper(path{i}(1)));
        lgc(i)=drive>=65 && drive<=90 && strcmp(path{i}(2:3),':\');
    end
else % unix, linux, macosx
    for i=1:numel(path)
        if(isempty(path{i})); continue; end
        lgc(i)=strcmp(path{i}(1),'/');
    end
end

end
