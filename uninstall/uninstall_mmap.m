function [ok]=uninstall_mmap()
%UNINSTALL_MMAP    Uninstalls the currently installed M_Map
%
%    Usage:    ok=uninstall_mmap
%
%    Description:
%     OK=UNINSTALL_MMAP removes the directory associated with the M_Map
%     toolbox from the path and tries to save the edited path.  OK is TRUE
%     if the uninstall succeeded.
%
%    Notes:
%     - Uses the location of function m_coast to detect the M_Map toolbox.
%
%    Examples:
%     % Update M_Map:
%     uninstall_mmap & webinstall_mmap
%
%    See also: WEBINSTALL_MMAP, UNINSTALL_GSHHS, WEBINSTALL_GSHHS,
%              UNINSTALL_NJTBX, WEBINSTALL_NJTBX, UNINSTALL_EXPORTFIG,
%              WEBINSTALL_EXPORTFIG, UNINSTALL_EXTRAS, WEBINSTALL_EXTRAS,
%              UNINSTALL_IRISWS, WEBINSTALL_IRISWS, UNINSTALL_TAUP,
%              WEBINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - handle not installed, flip logic from savepath,
%                        doc update
%        Apr. 25, 2012 - uninstall m_map_fixes directory too
%        Jan. 15, 2014 - updated See also list
%        Mar.  1, 2014 - updated See also list, savepath only called if
%                        needed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:25 GMT

% todo:

% does m_coast exist?
ok=true;
if(exist('m_coast','file'))
    path{1}=fileparts(which('m_coast')); % m_map directory
    path{2}=fileparts(which('m_gshhs')); % m_map_fixes directory
    rmpath(path{:});
    if(is_on_static_path(path{:}))
        ok=~savepath;
    end
else
    % not found, so toolbox not installed...
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

