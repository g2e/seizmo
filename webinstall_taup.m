function [ok]=webinstall_taup(mypath)
%WEBINSTALL_TAUP    Install TauP components
%
%    Usage:    ok=webinstall_taup
%              ok=webinstall_taup(path)
%
%    Description:
%     OK=WEBINSTALL_TAUP downloads & installs the jar files needed by
%     MatTauP into the 'mattaup/lib' directory.  The download is not large
%     at about 1 megabyte but do make sure you have a good connection or
%     this operation will take a little while and you cannot do anything
%     else at the Matlab/Octave prompt while waiting for the files to
%     download.
%
%     OK=WEBINSTALL_TAUP(PATH) installs the files in the directory given
%     by PATH.
%
%    Notes:
%     - Also adds the MatTauP-*.jar to the classpath.
%
%    Examples:
%     % Reinstall the MatTauP & TauP jar files:
%     uninstall_taup & webinstall_taup
%
%    See also: UNINSTALL_IRISWS, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_MMAP, UNINSTALL_MMAP, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXPORTFIG, WEBINSTALL_NJTBX, UNINSTALL_NJTBX,
%              WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, INSTALL_SEIZMO, UNINSTALL_SEIZMO

%     Version History:
%        Feb. 20, 2014 - initial version
%        Feb. 25, 2014 - bugfix: assign to output
%        Feb. 28, 2014 - install on dynamic & static (if possible)
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
    mypath=[fileparts(mfilename('fullpath')) fs 'mattaup' fs 'lib'];
    if(~exist(mypath,'dir')); mkdir(mypath); end
else
    if(~isabspath(mypath)); mypath=[pwd fs mypath]; end
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_taup:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt TauP install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    taup='TauP-2.1.1.jar';       % ~700kb
    seis='seisFile-1.5.1.jar';   % ~200kb
    mattaup='MatTauP-2.1.1.jar'; %  ~10kb (comes with SEIZMO)
    
    % grab files (either locally or remotely)
    url0='http://www.seis.sc.edu/software/maven2/edu/sc/seis/TauP/2.1.1/';
    url1=['http://www.seis.sc.edu/software/maven2/' ...
        'edu/sc/seis/seisFile/1.5.1/'];
    fprintf(' Getting %s\n',taup);
    if(exist(taup,'file'))
        if(~exist([mypath fs taup],'file'))
            copyfile(which(taup),'.');
        end
    else
        urlwrite([url0 taup],taup);
    end
    fprintf(' Getting %s\n',seis);
    if(exist(seis,'file'))
        if(~exist([mypath fs seis],'file'))
            copyfile(which(seis),'.');
        end
    else
        urlwrite([url1 seis],seis);
    end
    
    % check that java pkg is installed
    java_in_octave=true;
    if(exist('OCTAVE_VERSION','builtin')==5 && isempty(ver('java')))
        warning('seizmo:webinstall_taup:noJavaTbx',...
            'Java package missing from Octave!');
        java_in_octave=false;
    end
    
    % install jars to dynamic classpath
    taupjar=[mypath fs taup];
    seisjar=[mypath fs seis];
    mattaupjar=[mypath fs mattaup];
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),taup))))
        javaaddpath(taupjar);
    end
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),seis))))
        javaaddpath(seisjar);
    end
    if(java_in_octave && ...
            all(cellfun('isempty',strfind(javaclasspath('-all'),mattaup))))
        javaaddpath(mattaupjar);
    end
    
    % install jars to static classpath
    sjcp=which('classpath.txt');
    if(~isempty(sjcp))
        % read classpath.txt
        s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
        
        % detect offending classpath.txt lines
        injcp(1)=any(~cellfun('isempty',strfind(s2,taupjar)));
        injcp(2)=any(~cellfun('isempty',strfind(s2,seisjar)));
        injcp(3)=any(~cellfun('isempty',strfind(s2,mattaupjar)));
        
        % only add if not there already
        if(sum(injcp)<3)
            fid=fopen(sjcp,'a+');
            if(fid~=-1)
                fseek(fid,0,'eof');
                if(~injcp(1)); fprintf(fid,'%s\n',taupjar); end
                if(~injcp(2)); fprintf(fid,'%s\n',seisjar); end
                if(~injcp(3)); fprintf(fid,'%s\n',mattaupjar); end
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
