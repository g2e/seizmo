function [ok]=uninstall_gshhg()
%UNINSTALL_GSHHG    Uninstalls the currently installed GSHHG
%
%    Usage:    ok=uninstall_gshhg
%
%    Description:
%     OK=UNINSTALL_GSHHG removes the directory containing GSHHG binary
%     files from the path and tries to save the edited path.  OK is TRUE if
%     the uninstall succeeded.
%
%    Notes:
%     - Does not delete the directory or the downloaded zip file.
%     - Uses the location of binary file gshhs_c.b to detect the directory.
%
%    Examples:
%     % Update GSHHG:
%     uninstall_gshhg & webinstall_gshhg
%
%    See also: WEBINSTALL_GSHHG, UNINSTALL_MMAP, WEBINSTALL_MMAP,
%              UNINSTALL_NJTBX, WEBINSTALL_NJTBX, UNINSTALL_EXPORTFIG,
%              WEBINSTALL_EXPORTFIG, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - handle not installed, flip logic from savepath,
%                        doc update
%        Jan. 14, 2014 - update for gshhg
%        Jan. 15, 2014 - updated See also list
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 15, 2014 at 15:25 GMT

% todo:

% does nj_time exist?
if(exist('gshhs_c.b','file'))
    path=fileparts(which('gshhs_c.b')); % root directory
    rmpath(path);
    ok=~savepath;
else
    % not found, so not installed...
    ok=true;
    return;
end

end
