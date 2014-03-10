function [ok]=webinstall_njtbx(mypath)
%WEBINSTALL_NJTBX    Install njTBX components
%
%    Usage:    ok=webinstall_njtbx
%              ok=webinstall_njtbx(path)
%
%    Description:
%     OK=WEBINSTALL_NJTBX creates the directory 'njtbx' where this mfile is
%     located, moves into it and downloads & installs the njTBX toolbox
%     files.  The download is large (~30 megabytes) so make sure you have a
%     good connection.  Also be aware that you cannot do anything else at
%     the Matlab/Octave prompt while waiting for the files to download.
%
%     OK=WEBINSTALL_NJTBX(PATH) installs the njTBX toolbox in the directory
%     given by PATH.
%
%    Notes:
%
%    Examples:
%     % Update njTBX:
%     uninstall_njtbx & webinstall_njtbx
%
%    See also: UNINSTALL_NJTBX, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_MMAP, UNINSTALL_MMAP, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXPORTFIG, WEBINSTALL_IRISWS, UNINSTALL_IRISWS,
%              WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO, SZ_TOC_WW3

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - doc update, flip savepath logic, only use
%                        javaaddpath or edit classpath as needed
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Mar.  8, 2012 - minor code changes for clarity
%        Apr. 25, 2012 - use zipped version of the svn checkout of njtbx as
%                        the zipfile available via their website is broken
%        Jan. 14, 2014 - throw warning rather than error if problem
%        Jan. 15, 2014 - updated See also list
%        Jan. 27, 2014 - added isabspath for abs path fix to path option,
%                        added handling of octave without java
%        Feb. 20, 2014 - fixed warning id, update see also list
%        Feb. 27, 2014 - install on dynamic & static (if possible)
%        Mar.  2, 2014 - avoid warnings by search javaclasspath by filename
%        Mar. 10, 2014 - use java 1.5 for compatibility
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 10, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% directory separator
fs=filesep;

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=[fileparts(mfilename('fullpath')) fs 'njtbx'];
    if(~exist(mypath,'dir')); mkdir(mypath); end
else
    if(~isabspath(mypath)); mypath=[pwd fs mypath]; end
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_njtbx:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt njTBX install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    njtbx='njToolbox-2.0.zip';           % ~5mb
    toolsui='toolsUI-4.0.49.jar';        % ~20mb
    njtools='njTools-2.0.12_jre1.5.jar'; % ~5mb
    
    % grab files (either locally or remotely)
    % - remote location is from japan b/c sourceforge
    %   does not have direct links on the english site
    url0='http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/';
    url=['http://es.sourceforge.jp/frs/g_redir.php?m=jaist&f=%2F' ...
        'njtbx%2FnjTBX-downloads%2F'];
    fprintf(' Getting %s\n',njtbx);
    if(exist(njtbx,'file'))
        if(~exist([mypath fs njtbx],'file'))
            copyfile(which(njtbx),'.');
        end
    else
        urlwrite([url0 njtbx],njtbx);
    end
    fprintf(' Getting %s\n',toolsui);
    if(exist(toolsui,'file'))
        if(~exist([mypath fs toolsui],'file'))
            copyfile(which(toolsui),'.');
        end
    else
        urlwrite([url toolsui],toolsui);
    end
    fprintf(' Getting %s\n',njtools);
    if(exist(njtools,'file'))
        if(~exist([mypath fs njtools],'file'))
            copyfile(which(njtools),'.');
        end
    else
        urlwrite([url njtools],njtools);
    end
    
    % delete pre-existing directory if in Octave
    njtbxdir=[mypath fs njtbx(1:end-4)]; % strip .zip
    if(exist(njtbxdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',njtbxdir);
        y=rmdir(njtbxdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install njTBX
    unzip(njtbx);
    addpath(njtbxdir,...
        [njtbxdir fs 'examples'],...
        [njtbxdir fs 'njFunc'],...
        [njtbxdir fs 'njTBX-2.0'],...
        [njtbxdir fs 'njTBX-2.0' fs 'Utilities']);
    savepath;
    
    % check that java pkg is installed
    java_in_octave=true;
    if(exist('OCTAVE_VERSION','builtin')==5 && isempty(ver('java')))
        warning('seizmo:webinstall_njtbx:noJavaTbx',...
            ['Java package missing from Octave! ' ...
            'NJTBX cannot be installed!']);
        java_in_octave=false;
    end
    
    % install jars to dynamic classpath
    toolsuijar=[mypath fs toolsui];
    njtoolsjar=[mypath fs njtools];
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),toolsui))))
        javaaddpath(toolsuijar);
    end
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),njtools))))
        javaaddpath(njtoolsjar);
    end
    
    % install jars to static classpath
    sjcp=which('classpath.txt');
    if(~isempty(sjcp))
        % read classpath.txt
        s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
        
        % detect offending classpath.txt lines
        injcp(1)=any(~cellfun('isempty',strfind(s2,toolsuijar)));
        injcp(2)=any(~cellfun('isempty',strfind(s2,njtoolsjar)));
        
        % only add if not there already
        if(sum(injcp)<2)
            fid=fopen(sjcp,'a+');
            if(fid~=-1)
                fseek(fid,0,'eof');
                if(~injcp(1)); fprintf(fid,'%s\n',toolsuijar); end
                if(~injcp(2)); fprintf(fid,'%s\n',njtoolsjar); end
                fclose(fid);
            end
        end
    end
    
    % return
    cd(cwd);
    
    % all good?
    ok=java_in_octave;
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
