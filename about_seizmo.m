function []=about_seizmo()
%ABOUT_SEIZMO    Pop-Up display with SEIZMO Toolbox info
%
%    Usage:    about_seizmo
%
%    Description:
%     ABOUT_SEIZMO displays the version info of SEIZMO in a modal dialog
%     box.  This is basically just for installation fun.
%
%    Notes:
%
%    Examples:
%     % Provide a quick SEIZMO "splash" to the user:
%     about_seizmo
%
%    See also: INSTALL_SEIZMO, VER

%     Version History:
%        Dec. 30, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 30, 2010 at 15:25 GMT

% todo:
% - how to display hyperlinks?
% - truely awesome seizmo logo
% - truely awesome seizmo song
%   - kill song when popup gone
% - no button (ie no msgbox)

% retreive all necessary data
info=ver('seizmo');

% quick return if seizmo not installed
if(isempty(info)); return; end

% make popup msg
msg=sprintf([...
    'Package: ' info.Name '\n' ...
    'Version: ' info.Version '\n' ...
    'Release: ' info.Release '\n' ...
    'Date:    ' info.Date '\n']);
%    '<a href="matlab:help seizmo">Command Line Help</a>\n' ...
%    '<a href="http://epsc.wustl.edu/~ggeuler/codes/m/seizmo">SEIZMO Website</a>']);

% attempt icon retreival (otherwise create one)
try
    data=load(['seizmo_icon_' info.Release '.mat'],'icon');
catch
    % something spiffy here
    data.icon=zeros(64,64,3);
end

% the popup
msgbox(msg,info.Name,'custom',data.icon,hot(64),'modal');

end
