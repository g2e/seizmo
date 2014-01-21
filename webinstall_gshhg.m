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
%              WEBINSTALL_EXPORTFIG, UNINSTALL_SEIZMO, INSTALL_SEIZMO,
%              MMAP, M_MAP, M_GSHHS

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 15, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=fileparts(mfilename('fullpath'));
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
        if(~exist(fullfile(mypath,gshhg),'file'))
            copyfile(which(gshhg),'.');
        end
    else
        urlwrite(url,gshhg);
    end
    
    % delete pre-existing directory
    % - this is only for octave and might not be necessary anymore
    gshhgdir=fullfile(mypath,'gshhg');
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
    ok=~savepath;
    if(~ok)
        warning('seizmo:webinstall_gshhg:noWritePathdef',...
            'Cannot save path!');
    end
    
    % return
    cd(cwd);
catch
    le=lasterror;
    warning(le.identifier,le.message);
    ok=false;
    cd(cwd);
end

end
