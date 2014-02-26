function [ok]=uninstall_extras()
%UNINSTALL_EXTRAS    Uninstalls the extra SEIZMO components
%
%    Usage:    ok=uninstall_extras
%
%    Description:
%     OK=UNINSTALL_EXTRAS currently returns TRUE without doing anything as
%     the extra SEIZMO components are only downloaded and no adjustment is
%     made to the path.
%
%    Notes:
%
%    Examples:
%     % Force a redownload of extras:
%     webinstall_extras(true)
%
%    See also: WEBINSTALL_EXTRAS, UNINSTALL_NJTBX, WEBINSTALL_NJTBX,
%              UNINSTALL_MMAP, WEBINSTALL_MMAP, UNINSTALL_GSHHG,
%              WEBINSTALL_GSHHG, UNINSTALL_EXPORTFIG, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_IRISWS, WEBINSTALL_IRISWS, UNINSTALL_TAUP,
%              WEBINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 20, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 15:25 GMT

% todo:

% just return true
ok=true;

end
