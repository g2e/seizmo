function [ok]=uninstall_exportfig()
%UNINSTALL_EXPORTFIG    Uninstalls the currently installed export_fig
%
%    Usage:    ok=uninstall_exportfig
%
%    Description:
%     OK=UNINSTALL_EXPORTFIG removes the directory associated with the
%     export_fig toolbox from the path and tries to save the edited path.
%     OK is TRUE if the uninstall succeeded.
%
%    Notes:
%     - Uses the location of function 'export_fig' to detect the toolbox.
%
%    Examples:
%     % Update export_fig:
%     uninstall_exportfig & webinstall_exportfig
%
%    See also: WEBINSTALL_EXPORTFIG, UNINSTALL_GSHHS, WEBINSTALL_GSHHS,
%              UNINSTALL_NJTBX, WEBINSTALL_NJTBX, UNINSTALL_MMAP,
%              WEBINSTALL_MMAP

%     Version History:
%        Feb. 16, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2012 at 15:25 GMT

% todo:

% does export_fig exist?
if(exist('export_fig','file'))
    path=fileparts(which('export_fig')); % root directory
    rmpath(path);
    ok=~savepath;
else
    % not found, so toolbox not installed...
    ok=true;
    return;
end

end
