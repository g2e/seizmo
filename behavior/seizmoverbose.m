function [varargout]=seizmoverbose(varargin)
%SEIZMOVERBOSE    Turn verbose SEIZMO output on (TRUE) or off (FALSE)
%
%    Usage:    oldstate=seizmoverbose(newstate)
%
%    Description:
%     OLDSTATE=SEIZMOVERBOSE(NEWSTATE) toggles basic verbose messaging in
%     SEIZMO.  This may be a text-based progress bar and it may be some
%     simple description of what is happening.  Useful for debugging or
%     checking what SEIZMO is doing.  Input & output for SEIZMOVERBOSE are
%     logical scalars.
%
%    Notes:
%
%    Examples:
%     % Try add with verbose off:
%     seizmoverbose(false)
%     add(data,0);
%     seizmoverbose(true)
%
%    See also: SEIZMODEBUG, SEIZMORESET

%     Version History:
%        Dec.  2, 2009 - initial version
%        Jan. 29, 2010 - warn & fix if invalid value in SEIZMO
%        Jan. 28, 2012 - doc update
%        Mar.  1, 2012 - not verbose by default when parallel tbx in use
%        Mar. 15, 2012 - undo parallel detection (not needed)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 20:00 GMT

% todo:

% retrieve global
global SEIZMO

% get current value
if(nargout)
    try
        varargout{1}=SEIZMO.SEIZMOVERBOSE;
        if(~islogical(varargout{1}) || ~isscalar(varargout{1}))
            warning('seizmo:seizmoverbose:badState',...
                ['STATE of SEIZMOVERBOSE must be TRUE or FALSE!\n' ...
                'Using default verbosity (TRUE)!']);
            varargout{1}=true;
        end
    catch
        varargout{1}=true;
    end
end

% update value
if(nargin)
    if(~islogical(varargin{1}) || ~isscalar(varargin{1}))
        error('seizmo:seizmoverbose:badState',...
            'STATE of SEIZMOVERBOSE must be TRUE or FALSE!');
    end
    SEIZMO.SEIZMOVERBOSE=varargin{1};
end

end
