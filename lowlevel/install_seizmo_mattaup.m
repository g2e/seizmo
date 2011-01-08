function [ok]=install_seizmo_mattaup(mypath,varargin)
%INSTALL_SEIZMO_MATTAUP    Check & Install MatTauP for SEIZMO
%
%    Usage:    ok=install_seizmo_mattaup(mypath)
%              ok=install_seizmo_mattaup(mypath,type)
%
%    Description:
%     OK=INSTALL_SEIZMO_MATTAUP(MYPATH) installs the MatTauP files
%     associated with the current SEIZMO on the path and also adds the
%     matTaup.jar file to the global Matlab classpath.txt file.  The path
%     is saved to the pathdef.m file that the path was loaded from at
%     startup.
%
%     OK=INSTALL_SEIZMO_MATTAUP(MYPATH,TYPE) allows editing the pathdef.m
%     save preference.  See SAVEPATH_SEIZMO for details.  The default is
%     no input to SAVEPATH_SEIZMO.
%
%    Notes:
%     - SEIZMO's matTaup.jar is customized.  Use of matTaup.jar from
%       another source will likely cause issues!
%
%    Examples:
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MMAP, INSTALL_SEIZMO_WW3,
%              INSTALL_SEIZMO_CORE, INSTALL_SEIZMO_OPTIONAL

%     Version History:
%        Dec. 30, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 30, 2010 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default path input
if(nargin<1 || isempty(mypath)); mypath='.'; end

% check path
fs=filesep;
if(~exist(mypath,'dir'))
    error('seizmo:install_seizmo_mattaup:badPath',...
        ['SEIZMO directory (' mypath ') does not exist!']);
end

% check that classpath exists
ok=false;
sjcp=which('classpath.txt');
if(isempty(sjcp))
    warning('seizmo:install_seizmo_mattaup:noJavaClassPath',...
        ['SEIZMO will not be able to use MatTauP because there is no\n' ...
        'classpath.txt available!  If you are running Matlab, this\n' ...
        'is a problem!  For Octave users, this is normal and means\n' ...
        'that SEIZMO on Octave is missing features.']);
    return;
end

% check that classpath is synced
%s1=javaclasspath('-static');
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
%if(~isequal(s1,s2))
%    warning('seizmo:install_seizmo_mattaup:restartMatlab',...
%        ['The static javaclasspath & classpath.txt are out of sync!\n' ...
%        'Please restart Matlab before installing MatTauP!']);
%    return;
%end

% add seizmo mattaup mfiles to path
addpath([mypath fs 'mattaup']);
bad=savepath_seizmo(varargin{:});

% check if matTauP jar is installed, return if found
if(~isempty(cell2mat(strfind(s2,'matTaup.jar'))))
    warning('seizmo:install_seizmo_mattaup:previousInstall',...
        ['It appears MatTauP is already installed in classpath.txt!\n' ...
        'SEIZMO uses a modified version of the MatTaup.jar!\n' ...
        'Your installed version of MatTauP may cause SEIZMO to fail!']);
    return;
end

% this is the tough part (appending classpath.txt)
myjarpath=[mypath fs 'mattaup' fs 'lib' fs 'matTaup.jar'];
fid=fopen(sjcp,'a+');
if(fid<0)
    warning('seizmo:install_seizmo_mattaup:failedToOpen',...
        ['You must have Root/Administrator privileges to edit\n' ...
        'the classpath.txt file.  To fully install MatTauP\n' ...
        'you need to add the following line:\n' ...
        strrep(myjarpath,'\','\\') '\n' ...
        'to your Matlab''s classpath.txt located here:\n' ...
        strrep(sjcp,'\','\\') '\n\n' ...
        'This may be done by contacting your System Admin if you\n'
        'do not have Root/Administrator privileges.']);
else
    fseek(fid,0,'eof');
    fprintf(fid,'%s\n',myjarpath);
    fclose(fid);
    if(~bad); ok=true; end
end

end
