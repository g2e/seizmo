function [ok]=webinstall_gshhg(mypath)
%WEBINSTALL_GSHHG    Install GSHHG components
%
%    Usage:    ok=webinstall_gshhg
%              ok=webinstall_gshhg(path)
%
%    Description:
%     OK=WEBINSTALL_GSHHG downloads the GSHHG zip file to the directory
%     where this mfile is located, extracts its contents to the 'gshhg'
%     directory and adds/saves it to the path.  THE DOWNLOAD IS LARGE:
%     OVER 100 MEGABYTES & CAN TAKE HOURS ON A SLOW CONNECTION!  THERE IS
%     NO RESUME: IF IT FAILS, YOU WILL LOSE EVERYTHING YOU DOWNLOADED!  See
%     the notes below on how to install the GSHHG binaries manually.  Also
%     be aware that you cannot do anything else at the Matlab/Octave prompt
%     while waiting for the file to download.
%
%     OK=WEBINSTALL_GSHHG(PATH) install the GSHHG binary files under
%     the directory given by PATH.
%
%    Notes:
%     - Get the GSHHG binary files by downloading the archive from here:
%        ftp://ftp.soest.hawaii.edu/gshhg/gshhg-bin-2.2.4.zip
%       Extract the binaries from the archive using your unzip utility of
%       choice.
%     - If you have the GSHHG binary files extracted, just add the
%       directory containing them to your Matlab or Octave path using:
%        addpath directory/of/gshhg/
%       where the directory/of/gshhg needs to be changed to where the
%       GSHHG binary files are located on your computer.
%
%    Examples:
%     % GSHHG is big: this download can take _hours_.  This function is
%     % only useful if you have a fast connection!
%
%    See also: UNINSTALL_GSHHG, WEBINSTALL_MMAP, UNINSTALL_MMAP,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX, UNINSTALL_EXPORTFIG,
%              WEBINSTALL_EXPORTFIG, WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS,
%              WEBINSTALL_IRISWS, UNINSTALL_IRISWS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO, MMAP,
%              M_MAP, M_GSHHS

%     Version History:
%        Feb.  5, 2011 - initial version
%        Feb. 14, 2012 - renamed from seizmo_gshhs_webinstall to
%                        webinstall_gshhs, updated to use gshhs2 files, doc
%                        update, no savpath_seizmo, installs to location of
%                        this file under gshhg directory
%        Feb. 15, 2012 - doc update, flip savepath logic
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Apr.  8, 2013 - url moved (not updating to gshhg atm)
%        Jan. 14, 2014 - url moved, updated/renamed to latest gshhg, throw
%                        warning if problem encountered
%        Jan. 15, 2014 - updated See also list
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
    error('seizmo:webinstall_gshhg:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt GSHHG install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current version
    gshhg='gshhg-bin-2.2.4.zip';
    
    % grab file
    url='ftp://ftp.soest.hawaii.edu/gshhg/';
    url=[url gshhg];
    disp([' Getting ' gshhg]);
    if(exist(gshhg,'file'))
        if(~exist([mypath fs gshhg],'file'))
            copyfile(which(gshhg),'.');
        end
    else
        urlwrite(url,gshhg);
    end
    
    % delete pre-existing directory
    % - this is only for octave and might not be necessary anymore
    gshhgdir=[mypath fs 'gshhg'];
    if(exist(gshhgdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',gshhgdir);
        y=rmdir(gshhgdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install GSHHG
    unzip(gshhg,gshhgdir);
    addpath(gshhgdir);
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
