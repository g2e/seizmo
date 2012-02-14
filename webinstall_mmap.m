function [ok]=webinstall_mmap(mypath)
%WEBINSTALL_MMAP    Install m_map components
%
%    Usage:    ok=webinstall_mmap
%              ok=webinstall_mmap(path)
%
%    Description:
%     OK=WEBINSTALL_MMAP downloads the m_map zip file to the directory
%     where this mfile is located, extracts its contents to the 'm_map'
%     directory and adds/saves it to the path.  The download is only 0.6
%     megabytes so major worries.  Just be aware that you cannot do
%     anything else at the Matlab/Octave prompt while waiting for the file
%     to download.
%
%     OK=WEBINSTALL_MMAP(PATH) installs the m_map toolbox under the
%     directory given by PATH.
%
%    Notes:
%
%    Examples:
%     % Update m_map:
%     uninstall_mmap & webinstall_mmap
%
%    See also: UNINSTALL_MMAP, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_NJTBX, UNINSTALL_NJTBX

%     Version History:
%        Feb. 14, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2012 at 15:25 GMT

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

% attempt mmap install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current version
    mmap='m_map1.4.zip';
    
    % grab file
    url='http://www.eos.ubc.ca/%7Erich/';
    if(exist(mmap,'file'))
        if(~exist(fullfile(mypath,mmap),'file'))
            copyfile(which(mmap),'.');
        end
    else
        urlwrite([url mmap],mmap);
    end
    
    % unpack and install mmap
    unzip(mmap);
    addpath('m_map');
    ok=savepath;
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
