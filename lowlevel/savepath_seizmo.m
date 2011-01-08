function [bad]=savepath_seizmo(varargin)
%SAVEPATH_SEIZMO    Savepath with simple global/local path option
%
%    Usage:    savepath_seizmo
%              savepath_seizmo(path)
%              savepath_seizmo('global'|'local'|'either')
%
%    Description:
%     SAVEPATH_SEIZMO or SAVEPATH_SEIZMO(PATH) behaves exactly like the
%     function SAVEPATH.  Only 3 keywords are excepted below.
%
%     SAVEPATH_SEIZMO('GLOBAL'|'LOCAL'|'EITHER') will attempt to save to
%     the global 'pathdef.m' file or a local one in the current directory.
%     The 'either' option will first try to save to the global file then
%     the local one if that fails.
%
%    Notes:
%
%    Examples:
%
%    See also: SAVEPATH, INSTALL_SEIZMO, UNINSTALL_SEIZMO

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  1, 2011 - improve handling of no arg case
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  1, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% handle no input specially
if(~nargin)
    bad=savepath;
    return;
end

% check type
valid={'global' 'local' 'either'};
if(~ischar(type))
    error('seizmo:uninstall_seizmo:badInput',...
        ['Saved PATH location preference must be of the following:\n' ...
        sprintf('''%s'' ',valid{:}) '\n' ...
        'or a valid path!']);
end

% save path
switch lower(type)
    case 'global'
        % attempt root level save
        bad=savepath([matlabroot fs 'toolbox' fs 'local' fs 'pathdef.m']);
    case 'local'
        % attempt local level save
        bad=savepath('pathdef.m');
    case 'either'
        % attempt root level save
        bad=savepath([matlabroot fs 'toolbox' fs 'local' fs 'pathdef.m']);
        
        if(bad)
            % attempt local level save
            bad=savepath('pathdef.m');
        end
    otherwise
        bad=savepath(type);
end

end
