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
%     % Update M_Map:
%     uninstall_mmap & webinstall_mmap
%
%    See also: UNINSTALL_MMAP, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - add m_map_fixes, doc update, flip savepath logic
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Apr. 25, 2012 - copy m_map/private to m_map_fixes/private
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 25, 2012 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default path to seizmo directory
if(nargin<1 || isempty(mypath))
    mypath=fileparts(mfilename('fullpath'));
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
        if(~exist(fullfile(mypath,mmap),'file'))
            copyfile(which(mmap),'.');
        end
    else
        urlwrite([url mmap],mmap);
    end
    
    % delete pre-existing directory if in Octave
    mmapdir=fullfile(mypath,'m_map');
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
    if(exist(fullfile(mypath,'m_map_fixes'),'dir'))
        addpath(fullfile(mypath,'m_map_fixes'));
        copyfile(fullfile(mmapdir,'private'),...
            fullfile(mypath,'m_map_fixes','private'));
    end
    ok=~savepath;
    if(~ok)
        warning('seizmo:webinstall_mmap:noWritePathdef',...
            'Cannot save path!');
    end
    
    % return
    cd(cwd);
catch
    ok=false;
    cd(cwd);
end

end
