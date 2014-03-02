function [ok]=webinstall_irisws(mypath)
%WEBINSTALL_IRISWS    Install IRIS Web Services components
%
%    Usage:    ok=webinstall_irisws
%              ok=webinstall_irisws(path)
%
%    Description:
%     OK=WEBINSTALL_IRISWS downloads & installs the files distributed by
%     IRIS for coupling Matlab/Octave with IRIS web services into the 'ws'
%     directory.  The download is not large at <1 megabyte but do make sure
%     you have a good connection or this operation will take a little while
%     and you cannot do anything else at the Matlab/Octave prompt while
%     waiting for the files to download.
%
%     OK=WEBINSTALL_IRISWS(PATH) installs the files in the directory given
%     by PATH.
%
%    Notes:
%
%    Examples:
%     % Reinstall IRIS web services:
%     uninstall_irisws & webinstall_irisws
%
%    See also: UNINSTALL_IRISWS, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_MMAP, UNINSTALL_MMAP, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXPORTFIG, WEBINSTALL_NJTBX, UNINSTALL_NJTBX,
%              WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, INSTALL_SEIZMO, UNINSTALL_SEIZMO

%     Version History:
%        Feb. 20, 2014 - initial version
%        Feb. 25, 2014 - use 'ws' directory rather than 'irisws',
%                        bugfix: assign to output
%        Feb. 27, 2014 - install on dynamic & static (if possible)
%        Mar.  2, 2014 - avoid warnings by search javaclasspath by filename
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  2, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% directory separator
fs=filesep;

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=[fileparts(mfilename('fullpath')) fs 'ws'];
    if(~exist(mypath,'dir')); mkdir(mypath); end
else
    if(~isabspath(mypath)); mypath=[pwd fs mypath]; end
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_irisws:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt IRISWS install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    irisws='IRIS-WS-2.0.6.jar'; % ~500kb
    irisfetch='irisFetch.m';    % ~100kb
    
    % grab files (either locally or remotely)
    url0='http://www.iris.edu/files/IRIS-WS/2/2.0.6/';
    url1='http://www.iris.edu/files/irisFetch.m/2-0-6/';
    fprintf(' Getting %s\n',irisws);
    if(exist(irisws,'file'))
        if(~exist([mypath fs irisws],'file'))
            copyfile(which(irisws),'.');
        end
    else
        urlwrite([url0 irisws],irisws);
    end
    fprintf(' Getting %s\n',irisfetch);
    if(exist(irisfetch,'file'))
        if(~exist([mypath fs irisfetch],'file'))
            copyfile(which(irisfetch),'.');
        end
    else
        urlwrite([url1 irisfetch],irisfetch);
    end
    
    % check that java pkg is installed
    java_in_octave=true;
    if(exist('OCTAVE_VERSION','builtin')==5 && isempty(ver('java')))
        warning('seizmo:webinstall_irisws:noJavaTbx',...
            ['Java package missing from Octave! ' ...
            'IRISWS cannot be installed']);
        java_in_octave=false;
    end
    
    % install jars to dynamic classpath if not already on dynamic/static
    % - Avoid having jars of the same name with different paths as this
    %   will generate warnings from at least Matlab.
    iriswsjar=[mypath fs irisws];
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),irisws))))
        javaaddpath(iriswsjar);
    end
    
    % install jars to static classpath
    sjcp=which('classpath.txt');
    if(~isempty(sjcp))
        % read classpath.txt
        s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
        
        % detect offending classpath.txt lines
        injcp(1)=any(~cellfun('isempty',strfind(s2,iriswsjar)));
        
        % only add if not there already
        if(sum(injcp)<1)
            fid=fopen(sjcp,'a+');
            if(fid~=-1)
                fseek(fid,0,'eof');
                if(~injcp(1)); fprintf(fid,'%s\n',iriswsjar); end
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
