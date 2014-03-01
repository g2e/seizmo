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
%              WEBINSTALL_EXPORTFIG, UNINSTALL_EXTRAS, WEBINSTALL_EXTRAS,
%              UNINSTALL_IRISWS, WEBINSTALL_IRISWS, UNINSTALL_TAUP,
%              WEBINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - handle not installed, flip logic from savepath,
%                        doc update
%        Jan. 14, 2014 - update for gshhg
%        Jan. 15, 2014 - updated See also list
%        Mar.  1, 2014 - updated See also list, savepath only called if
%                        needed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:25 GMT

% todo:

% does nj_time exist?
ok=true;
if(exist('gshhs_c.b','file'))
    path=fileparts(which('gshhs_c.b')); % root directory
    rmpath(path);
    if(is_on_static_path(path))
        ok=~savepath;
    end
else
    % not found, so not installed...
    return;
end

end


function [lgc]=is_on_static_path(varargin)
% find pathdef.m
spd=which('pathdef.m');

% read pathdef.m
s=textread(spd,'%s','delimiter','\n','whitespace','');

% detect offending pathdef.m lines
for i=1:nargin
    lgc=any(~cellfun('isempty',strfind(s,varargin{i})));
    if(lgc); return; end
end
end

