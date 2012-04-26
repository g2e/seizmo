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
%    See also: UNINSTALL_NJTBX, UNINSTALL_GSHHS, WEBINSTALL_GSHHS,
%              UNINSTALL_MMAP, WEBINSTALL_MMAP

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - doc update, flip savepath logic, only use
%                        javaaddpath or edit classpath as needed
%        Feb. 16, 2012 - workaround quietly stalled unzip in octave
%        Mar.  8, 2012 - minor code changes for clarity
%        Apr. 25, 2012 - use zipped version of the svn checkout of njtbx as
%                        the zipfile available via their website is broken
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 25, 2012 at 15:25 GMT

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

% attempt njTBX install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    njtbx='njToolbox-2.0.zip';           % ~5mb
    toolsui='toolsUI-4.0.49.jar';        % ~20mb
    njtools='njTools-2.0.12_jre1.6.jar'; % ~5mb
    
    % grab files (either locally or remotely)
    % - remote location is from japan b/c sourceforge
    %   does not have direct links on the english site
    url0='http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/';
    url=['http://es.sourceforge.jp/frs/g_redir.php?m=jaist&f=%2F' ...
        'njtbx%2FnjTBX-downloads%2F'];
    fprintf(' Getting %s\n',njtbx);
    if(exist(njtbx,'file'))
        if(~exist(fullfile(mypath,njtbx),'file'))
            copyfile(which(njtbx),'.');
        end
    else
        urlwrite([url0 njtbx],njtbx);
    end
    fprintf(' Getting %s\n',toolsui);
    if(exist(toolsui,'file'))
        if(~exist(fullfile(mypath,toolsui),'file'))
            copyfile(which(toolsui),'.');
        end
    else
        urlwrite([url toolsui],toolsui);
    end
    fprintf(' Getting %s\n',njtools);
    if(exist(njtools,'file'))
        if(~exist(fullfile(mypath,njtools),'file'))
            copyfile(which(njtools),'.');
        end
    else
        urlwrite([url njtools],njtools);
    end
    
    % delete pre-existing directory if in Octave
    njtbxdir=fullfile(mypath,njtbx(1:end-4)); % strip .zip
    if(exist(njtbxdir,'dir') && exist('OCTAVE_VERSION','builtin')==5)
        fprintf('Output directory exists: %s\n',njtbxdir);
        y=rmdir(njtbxdir,'s');
        if(~y)
            disp('Replace All or None of the files? A/N?');
        end
    end
    
    % unpack and install njTBX
    unzip(njtbx);
    fs=filesep;
    addpath(njtbxdir,...
        [njtbxdir fs 'examples'],...
        [njtbxdir fs 'njFunc'],...
        [njtbxdir fs 'njTBX-2.0'],...
        [njtbxdir fs 'njTBX-2.0' fs 'Utilities']);
    ok=~savepath;
    if(~ok)
        warning('seizmo:webinstall_njtbx:noWritePathdef',...
            'Cannot save path!');
    end
    
    % install jars to classpath
    toolsuijar=fullfile(mypath,toolsui);
    njtoolsjar=fullfile(mypath,njtools);
    sjcp=which('classpath.txt');
    if(isempty(sjcp))
        %warning('seizmo:webinstall_njtbx:noJavaClassPath',...
        %    'Octave has no classpath.txt to save .jar files!');
        
        % no classpath.txt so add to dynamic path
        if(~ismember(toolsuijar,javaclasspath))
            javaaddpath(toolsuijar);
        end
        if(~ismember(njtoolsjar,javaclasspath))
            javaaddpath(njtoolsjar);
        end
    else
        % read classpath.txt
        s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
        
        % detect offending classpath.txt lines
        injcp(1)=any(~cellfun('isempty',strfind(s2,toolsuijar)));
        injcp(2)=any(~cellfun('isempty',strfind(s2,njtoolsjar)));
        
        % only add if not there already
        if(sum(injcp)<2)
            fid=fopen(sjcp,'a+');
            if(fid<0)
                warning('seizmo:webinstall_njtbx:noWriteClasspath',...
                    ['Cannot edit classpath.txt! Adding njTBX jars ' ...
                    'to dynamic java class path!']);
                if(~ismember(toolsuijar,javaclasspath))
                    javaaddpath(toolsuijar);
                end
                if(~ismember(njtoolsjar,javaclasspath))
                    javaaddpath(njtoolsjar);
                end
            else
                fseek(fid,0,'eof');
                if(~injcp(1)); fprintf(fid,'%s\n',toolsuijar); end
                if(~injcp(2)); fprintf(fid,'%s\n',njtoolsjar); end
                fclose(fid);
            end
        end
    end
    
    % return
    cd(cwd);
catch
    error(lasterror)
    ok=false;
    cd(cwd);
end

end
