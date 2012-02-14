function [ok]=webinstall_njtbx(mypath)
%WEBINSTALL_NJTBX    Install njtbx components
%
%    Usage:    ok=webinstall_njtbx
%              ok=webinstall_njtbx(path)
%
%    Description:
%     OK=WEBINSTALL_NJTBX creates the directory 'njtbx' where this mfile is
%     located, moves into it and downloads & installs the njtbx toolbox
%     files.  The download is large (~30 megabytes) so make sure you have a
%     good connection.  Also be aware that you cannot do anything else at
%     the Matlab/Octave prompt while waiting for the files to download.
%
%     OK=WEBINSTALL_NJTBX(PATH) installs the njtbx toolbox in the directory
%     given by PATH.
%
%    Notes:
%
%    Examples:
%     % Update njtbx:
%     uninstall_njtbx & webinstall_njtbx
%
%    See also: UNINSTALL_NJTBX, UNINSTALL_GSHHS, WEBINSTALL_GSHHS,
%              UNINSTALL_MMAP, WEBINSTALL_MMAP

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
    mypath=fullfile(fileparts(mfilename('fullpath')),'njtbx','');
    if(~exist(mypath,'dir')); mkdir(mypath); end
end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_njtbx:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% check that classpath exists (Octave fails here)
sjcp=which('classpath.txt');
if(isempty(sjcp))
    warning('seizmo:webinstall_njtbx:noJavaClassPath',...
        'Octave has no classpath.txt to save .jar files!');
end

% attempt njtbx install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    njtbx='matlab-njTbx-2.0.05.zip';     % ~5mb
    toolsui='toolsUI-4.0.49.jar';        % ~20mb
    njtools='njTools-2.0.12_jre1.6.jar'; % ~5mb
    
    % grab files (either locally or remotely)
    % - remote location is from japan b/c sourceforge
    %   does not have direct links on the english site
    url=['http://es.sourceforge.jp/frs/g_redir.php?m=jaist&f=%2F' ...
        'njtbx%2FnjTBX-downloads%2F'];
    disp(['Getting ' njtbx]);
    if(exist(njtbx,'file'))
        if(~exist(fullfile(mypath,njtbx),'file'))
            copyfile(which(njtbx),'.');
        end
    else
        urlwrite([url njtbx],njtbx);
    end
    disp(['Getting ' toolsui]);
    if(exist(toolsui,'file'))
        if(~exist(fullfile(mypath,toolsui),'file'))
            copyfile(which(toolsui),'.');
        end
    else
        urlwrite([url toolsui],toolsui);
    end
    disp(['Getting ' njtools]);
    if(exist(njtools,'file'))
        if(~exist(fullfile(mypath,njtools),'file'))
            copyfile(which(njtools),'.');
        end
    else
        urlwrite([url njtools],njtools);
    end
    
    % unpack and install njtbx
    unzip(njtbx);
    addpath(fullfile(njtbx,'njTBX-2.0','Utilities'));
    addpath(fullfile(njtbx,'njTBX-2.0'));
    addpath(fullfile(njtbx,'njFunc'));
    addpath(fullfile(njtbx,'examples'));
    addpath(njtbx);
    ok=savepath;
    if(~ok)
        warning('seizmo:webinstall_njtbx:noWritePathdef',...
            'Cannot save path!');
    end
    
    % install jars to classpath
    if(isempty(sjcp))
        % no classpath.txt so add to dynamic path
        javaaddpath(toolsui);
        javaaddpath(njtools);
        ok=false;
    else
        fid=fopen(sjcp,'a+');
        if(fid<0)
            warning('seizmo:webinstall_njtbx:noWriteClasspath',...
                'Cannot edit classpath.txt!');
            ok=false;
        else
            fseek(fid,0,'eof');
            fprintf(fid,'%s\n',toolsui);
            fprintf(fid,'%s\n',njtools);
            fclose(fid);
            ok=ok & true;
        end
    end
    
    % return
    cd(cwd);
catch
    ok=false;
    cd(cwd);
end

end
