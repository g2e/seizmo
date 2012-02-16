function [ok]=uninstall_gshhs()
%UNINSTALL_GSHHS    Uninstalls the currently installed gshhs
%
%    Usage:    ok=uninstall_gshhs
%
%    Description:
%     OK=UNINSTALL_GSHHS removes the directory associated with the gshhs
%     binary files from the path and tries to save the edited path.  OK is
%     the savepath exit status.
%
%    Notes:
%     - Uses the location of binary file gshhs_c.b to detect the directory.
%
%    Examples:
%     % Update gshhs:
%     uninstall_gshhs & webinstall_gshhs
%
%    See also: WEBINSTALL_GSHHS, UNINSTALL_MMAP, WEBINSTALL_MMAP,
%              UNINSTALL_NJTBX, WEBINSTALL_NJTBX

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - handle not installed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2012 at 15:25 GMT

% todo:

% does nj_time exist?
if(exist('gshhs_c.b','file'))
    path=fileparts(which('gshhs_c.b')); % root directory
    rmpath(path);
    ok=savepath;
else
    % not found, so not installed...
    ok=true;
    return;
end

end
