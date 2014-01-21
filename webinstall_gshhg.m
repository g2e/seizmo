function [ok]=webinstall_gshhs(mypath)
%WEBINSTALL_GSHHS    Install GSHHS components
%
%    Usage:    ok=webinstall_gshhs
%              ok=webinstall_gshhs(path)
%
%    Description:
%     OK=WEBINSTALL_GSHHS downloads the GSHHS zip file to the directory
%     where this mfile is located, extracts its contents to the 'gshhs'
%     directory and adds/saves it to the path.  THE DOWNLOAD IS LARGE:
%     OVER 100 MEGABYTES & CAN TAKE HOURS ON A SLOW CONNECTION!  THERE IS
%     NO RESUME: IF IT FAILS, YOU WILL LOSE EVERYTHING YOU DOWNLOADED!  See
%     the notes below on how to install the GSHHS binaries manually.  Also
%     be aware that you cannot do anything else at the Matlab/Octave prompt
%     while waiting for the file to download.
%
%     OK=WEBINSTALL_GSHHS(PATH) install the GSHHS binary files under
%     the directory given by PATH.
%
%    Notes:
%     - Get the GSHHS binary files by downloading the archive from here:
%        ftp://ftp.soest.hawaii.edu/pwessel/gshhs/gshhs+wdbii_2.2.0.zip
%       Extract the binaries from the archive using your unzip utility of
%       choice.
%     - If you have the GSHHS binary files extracted, just add the
%       directory containing them to your Matlab or Octave path using:
%        addpath directory/of/gshhs/
%       where the directory/of/gshhs needs to be changed to where the
%       GSHHS binary files are located on your computer.
%
%    Examples:
%     % GSHHS is big: this download can take _hours_.  This function is
%     % only useful if you have a fast connection!
%
%    See also: UNINSTALL_GSHHS, WEBINSTALL_MMAP, UNINSTALL_MMAP,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX, M_MAP, M_GSHHS

%     Version History:
%        Feb.  5, 2011 - initial version
%        Feb. 14, 2012 - renamed from seizmo_gshhs_webinstall to
%                        webinstall_gshhs, updated to use gshhs2 files, doc
%                        update, no savpath_seizmo, installs to location of
%                        this file under gshhs directory
%        Feb. 15, 2012 - doc update, flip savepath logic
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Apr.  8, 2013 - url moved (not updating to gshhg atm)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  8, 2013 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=fileparts(mfilename('fullpath'));
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_gshhs:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt GSHHS install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current version
    gshhs='gshhs+wdbii_2.2.0.zip';
    
    % grab file
    url='ftp://ftp.soest.hawaii.edu/pwessel/gshhs/';
    url=[url gshhs];
    disp([' Getting ' gshhs]);
    if(exist(gshhs,'file'))
        if(~exist(fullfile(mypath,gshhs),'file'))
            copyfile(which(gshhs),'.');
        end
    else
        urlwrite(url,gshhs);
    end
    
    % delete pre-existing directory
    gshhsdir=fullfile(mypath,'gshhs');
    if(exist(gshhsdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',gshhsdir);
        y=rmdir(gshhsdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install GSHHS
    unzip(gshhs);
    addpath(gshhsdir);
    ok=~savepath;
    if(~ok)
        warning('seizmo:webinstall_gshhs:noWritePathdef',...
            'Cannot save path!');
    end
    
    % return
    cd(cwd);
catch
    ok=false;
    cd(cwd);
end

end
