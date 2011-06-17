function [varargout]=install_seizmo(varargin)
%INSTALL_SEIZMO    Installs the SEIZMO Toolbox in Matlab or Octave
%
%    Usage:    install_seizmo
%              install_seizmo(type)
%
%    Description:
%     INSTALL_SEIZMO installs the SEIZMO toolbox in Matlab or Octave.  This
%     mainly involves editing the Matlab/Octave path so you get access to
%     all the functions that make SEIZMO what it is.  There are a few
%     extras that require more work: using the MatTauP toolset requires
%     editing the system-wide classpath.txt (use 'which classpath.txt' in
%     Matlab to locate this -- sorry Octave does not yet support this), the
%     WaveWatch III toolset requires compiling a mex file, and the M_Map
%     toolbox needs the GSHHS coastline files downloaded and put on the
%     path too.  While the first 2 of these tasks are done here, the GSHHS
%     download is currently left to the user (see SEIZMO_GSHHS_WEBINSTALL).
%     An imformative message is given during the installation to aid you in
%     completing the task.  Also, for those interested this will uninstall
%     previous SEIZMO installs.  You can do an uninstallation yourself by
%     calling UNINSTALL_SEIZMO (its hidden in the 'lowlevel' directory).
%
%     INSTALL_SEIZMO(TYPE) allows editing the pathdef.m save preference.
%     See SAVEPATH_SEIZMO for details.  The default is no input to
%     SAVEPATH_SEIZMO.
%
%    Notes:
%     - INSTALL_SEIZMO must be called from *WITHIN* a running session of
%       Matlab or Octave.  This is NOT a shell script or the Windows
%       equivalent, so do NOT "run" install_seizmo.m from your OS shell.
%       Start up Matlab/Octave and in the command window type
%       "install_seizmo" without the quotes and press enter.  Read what it
%       says and with any luck you will be ready to go!
%
%     - Websites:
%        SEIZMO  - http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/
%        MatTauP - http://www.ess.washington.edu/SEIS/FMI/matTaup.htm
%        M_Map   - http://www.eos.ubc.ca/~rich/map.html
%
%    Examples:
%     % Amazingly, every step of installing SEIZMO can be done *WITHIN*
%     % Matlab/Octave.  All you need is an internet connection so you
%     % can grab the SEIZMO package, extract the contents and install:
%     file=gunzip('https://github.com/g2e/seizmo/tarball/master');
%     untar(file{:});
%     cd seizmo/
%     install_seizmo
%
%     % and if you have the time:
%     seizmo_gshhs_webinstall
%
%    See also: ABOUT_SEIZMO, SEIZMO, UNINSTALL_SEIZMO, SAVEPATH_SEIZMO,
%              SEIZMO_GSHHS_WEBINSTALL

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  1, 2011 - added msgs, detects more issues
%        Jan.  4, 2011 - improved docs (examples needs more work)
%        Jan. 14, 2011 - example improved to get current release
%        Jan. 19, 2011 - updated to fixed examples
%        Apr.  6, 2011 - include verLessThan for pre-7.4 matlab
%        June 16, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% check application & version
disp('##################################################################');
disp('################## STARTING SEIZMO INSTALLATION ##################');
disp('##################################################################');
disp('Checking what application we are installing SEIZMO in...')
[application,version]=getapplication;
disp(['Application:  ' application]);
disp(['Version    :  ' version]);
disp(' ');
switch lower(application)
    case 'matlab'
        if(verLessThan('matlab','7.1'))
            warning('seizmo:install_seizmo:versionBad',...
                ['Matlab version too old for SEIZMO!\n' ...
                'Full function of toolbox is unlikely.']);
        end
        if(~license('checkout','signal_toolbox'))
            warning('seizmo:install_seizmo:noSigProcTbx',...
                ['Your Matlab does not have the Signal Processing\n' ...
                'Toolbox installed!  SEIZMO needs the Signal\n' ...
                'Processing Toolbox to be fully functional!']);
        elseif(~license('checkout','statistics_toolbox'))
            warning('seizmo:install_seizmo:noStatsTbx',...
                ['Your Matlab does not have the Statistics\n' ...
                'Toolbox installed!  SEIZMO needs the Statistics\n' ...
                'Toolbox to be fully functional!']);
        end
    case 'octave'
        warning('seizmo:install_seizmo:octaveIssues',...
            ['Octave compatibility for SEIZMO is a work in progress.\n' ...
            'Please report issues as found.  Java issues (MatTauP)\n' ...
            'and mex compilation problems are expected.']);
    otherwise
        warning('seizmo:install_seizmo:noClueWhatIsRunning',...
            'Unsure if SEIZMO will install on UNKNOWN application!');
end

% where am i?
me=mfilename('fullpath');
mypath=fileparts(me);
cwd=pwd;
disp(' ');
disp(['SEIZMO install path:  ' mypath]);
disp(' ');

% what is my version
disp('Looking for previous SEIZMO installs...');
info=ver('seizmo');

% remove old seizmo installations
ok=true(6,1);
while(~isempty(info))
    disp('Found previous SEIZMO installation.  Uninstalling!');
    ok(1)=uninstall_seizmo();
    if(~ok(1))
        warning('seizmo:install_seizmo:failedUninstall',...
            ['Uninstalling previous SEIZMO failed.  Please\n' ...
            'resolve this and then attempt INSTALL_SEIZMO again.\n']);
        if(nargout); varargout{1}=false; end
        return;
    end
    info=ver('seizmo');
end

% install new seizmo
if(ok(1))
    cd([mypath filesep 'lowlevel']);
    disp('Installing SEIZMO-MatTauP components...');
    ok(2)=install_seizmo_mattaup(mypath,varargin{:});
    if(~ok(2))
        warning('seizmo:install_seizmo:failedInstall',...
            'Failed to install SEIZMO-MatTauP components!');
    end
    disp('Installing SEIZMO-M_Map components...');
    ok(3)=install_seizmo_mmap(mypath,varargin{:});
    if(~ok(3))
        warning('seizmo:install_seizmo:failedInstall',...
            'Failed to install SEIZMO-M_Map components!');
    end
    disp('Installing SEIZMO-WaveWatch III components...');
    ok(4)=install_seizmo_ww3(mypath,varargin{:});
    if(~ok(4))
        warning('seizmo:install_seizmo:failedInstall',...
            'Failed to install SEIZMO-WaveWatch III components!');
    end
    disp('Installing Optional SEIZMO components...');
    ok(5)=install_seizmo_optional(mypath,varargin{:});
    if(~ok(5))
        warning('seizmo:install_seizmo:failedInstall',...
            'Failed to install optional SEIZMO components!');
    end
    disp('Installing Core SEIZMO components...');
    ok(6)=install_seizmo_core(mypath,varargin{:});
    if(~ok(6))
        warning('seizmo:install_seizmo:failedInstall',...
            'Failed to install core SEIZMO components!');
    end
    cd(cwd);
end

disp('##################################################################');
disp('################## FINISHED SEIZMO INSTALLATION ##################');
disp('##################################################################');
fprintf('\n\n\n');

% display helpful info
about_seizmo;
help seizmo;

% output
if(nargout); varargout{1}=sum(ok)==numel(ok); end

end

function [application,version]=getapplication()
%GETAPPLICATION    Returns application running this script and its version
%
%    Usage:    [application,version]=getapplication()
%
%    Description: [APPLICATION,VERSION]=GETAPPLICATION() will determine and
%     return the name and version of the application running this script
%     (obviously only if the application can run this script in the first
%     place).  Both APPLICATION and VERSION are strings.
%
%    Notes:
%     - returns 'UNKNOWN' if it cannot figure out the application
%
%    Examples:
%     Matlab and Octave still behave quite differently for a number of
%     different functions so it is best in some cases to use different
%     function calls depending on which we are running:
%      [app,ver]=getapplication;
%      if(strcmp(app,'MATLAB'))
%        % do something via matlab routines
%      else
%        % do something via octave routines
%      end
%
%    See also: NATIVEBYTEORDER, VER

%     Version History:
%        Nov. 13, 2008 - initial version
%        Mar.  3, 2009 - minor doc cleaning
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 19:55 GMT

% todo:

% checking for Matlab will throw an error in Octave
try
    % first check if we are in Matlab
    a=ver('matlab');
    
    % we are in Matlab
    application=a.Name;
    version=a.Version;
    return;
catch
    % check if we are in Octave
    if(exist('OCTAVE_VERSION','builtin')==5)
        application='OCTAVE';
        version=OCTAVE_VERSION;
        return;
    % ok I have no clue what is running
    else
        application='UNKNOWN';
        version='UNKNOWN';
        return;
    end
end

end

function result = verLessThan(toolboxstr, verstr)
%verLessThan Compare version of toolbox to specified version string.
%   verLessThan(TOOLBOX_DIR, VERSION) returns true if the version of
%   the toolbox specified by the string TOOLBOX_DIR is older than the
%   version specified by the string VERSION, and false otherwise. 
%   VERSION must be a string in the form 'major[.minor[.revision]]', 
%   such as '7', '7.1', or '7.0.1'. If TOOLBOX_DIR cannot be found
%   on MATLAB's search path, an error is generated.
%
%   Examples:
%       if verLessThan('images', '4.1')
%           error('Image Processing Toolbox 4.1 or higher is required.');
%       end
%
%       if verLessThan('matlab', '7.0.1')
%           % Put code to run under MATLAB older than MATLAB 7.0.1 here
%       else
%           % Put code to run under MATLAB 7.0.1 and newer here
%       end
%
%   See also MATLABPATH, VER.

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $  $Date: 2007/02/23 13:28:34 $
    
if nargin < 2
    errstr = 'Not enough input arguments.';
    if errorSupportsIdentifiers
        error('MATLAB:nargchk:notEnoughInputs', errstr)
    else
        error(errstr)
    end
end

if ~ischar(toolboxstr) | ~ischar(verstr)
    errstr = 'Inputs must be strings.';
    if errorSupportsIdentifiers
        error('MATLAB:verLessThan:invalidInput', errstr)
    else
        error(errstr)
    end
end

toolboxver = ver(toolboxstr);
if isempty(toolboxver)
    errformat = 'Toolbox ''%s'' not found.';
    if errorSupportsIdentifiers
        error('MATLAB:verLessThan:missingToolbox', errformat, toolboxstr)
    else
        error(sprintf(errformat, toolboxstr))
    end
end

toolboxParts = getParts(toolboxver(1).Version);
verParts = getParts(verstr);

result = (sign(toolboxParts - verParts) * [1; .1; .01]) < 0;
end

function parts = getParts(V)
    parts = sscanf(V, '%d.%d.%d')';
    if length(parts) < 3
       parts(3) = 0; % zero-fills to 3 elements
    end
end

function tf = errorSupportsIdentifiers
    % Determine, using code that runs on MATLAB 6.0 or later, if
    % error identifiers should be used when calling error().
    tf = 1;
    eval('lasterr('''','''');','tf = 0;');
end
