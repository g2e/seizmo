function [ok]=webinstall_exportfig(mypath)
%WEBINSTALL_EXPORTFIG    Install export_fig components
%
%    Usage:    ok=webinstall_exportfig
%              ok=webinstall_exportfig(path)
%
%    Description:
%     OK=WEBINSTALL_EXPORTFIG downloads the export_fig zip file to the
%     directory where this mfile is located, extracts its contents to the
%     'export_fig' directory and adds/saves it to the path.  The download
%     is only 0.1 megabytes.  Just be aware that you will not be able to do
%     anything else at the Matlab/Octave prompt while waiting for the file
%     to download.
%
%     OK=WEBINSTALL_EXPORTFIG(PATH) installs the export_fig toolbox under
%     the directory given by PATH.
%
%    Notes:
%
%    Examples:
%     % Update export_fig:
%     uninstall_exportfig & webinstall_exportfig
%
%    See also: UNINSTALL_EXPORTFIG, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX, WEBINSTALL_MMAP,
%              UNINSTALL_MMAP

%     Version History:
%        Feb. 16, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2012 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=fileparts(mfilename('fullpath'));
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_exportfig:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% attempt export_fig install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current version
    exportfig='export_fig.zip';
    
    % grab file
    url=['http://www.mathworks.com/matlabcentral/fileexchange/' ...
        '23629-exportfig?controller=file_infos&download=true'];
    disp([' Getting ' exportfig]);
    if(exist(exportfig,'file'))
        if(~exist(fullfile(mypath,exportfig),'file'))
            copyfile(which(exportfig),'.');
        end
    else
        urlwrite(url,exportfig);
    end
    
    % delete pre-existing directory if in Octave
    exportfigdir=fullfile(mypath,exportfig(1:end-4)); % strip .zip
    if(exist(exportfigdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',exportfigdir);
        y=rmdir(exportfigdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install export_fig
    unzip(exportfig,exportfig(1:end-4)); % extract to export_fig dir
    addpath(fullfile(mypath,'export_fig'));
    ok=~savepath;
    if(~ok)
        warning('seizmo:webinstall_exportfig:noWritePathdef',...
            'Cannot save path!');
    end
    
    % return
    cd(cwd);
catch
    ok=false;
    cd(cwd);
end

end
