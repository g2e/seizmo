function [varargout]=seizmodebug(varargin)
%SEIZMODEBUG    Turn SEIZMO debugging output on (TRUE) or off (FALSE)
%
%    Usage:    oldstate=seizmodebug(newstate)
%
%    Description:
%     OLDSTATE=SEIZMODEBUG(NEWSTATE) toggles debugging messages in SEIZMO.
%     This usually does nothing but may lead to some poorly-understood &
%     arcane output.  Useful for intense debugging.  Input & output for
%     SEIZMODEBUG are logical scalars.
%
%    Notes:
%
%    Examples:
%     % Try checkheader with debugging on:
%     seizmodebug(true)
%     checkheader(data);
%     seizmodebug(false)
%
%    See also: SEIZMOVERBOSE, SEIZMORESET

%     Version History:
%        Feb.  4, 2010 - initial version
%        Jan. 28, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 20:00 GMT

% todo:

% retrieve global
global SEIZMO

% get current value
if(nargout)
    try
        varargout{1}=SEIZMO.SEIZMODEBUG;
        if(~islogical(varargout{1}) || ~isscalar(varargout{1}))
            warning('seizmo:seizmodebug:badState',...
                ['STATE of SEIZMODEBUG must be TRUE or FALSE!\n' ...
                'Using default debug state (FALSE)!']);
            varargout{1}=false;
        end
    catch
        varargout{1}=false;
    end
end

% update value
if(nargin)
    if(~islogical(varargin{1}) || ~isscalar(varargin{1}))
            error('seizmo:seizmodebug:badState',...
                'STATE of SEIZMODEBUG must be TRUE or FALSE!');
    end
    SEIZMO.SEIZMODEBUG=varargin{1};
end

end
