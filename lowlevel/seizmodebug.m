function [varargout]=seizmodebug(varargin)
%SEIZMODEBUG    Turn SEIZMO debugging output on (TRUE) or off (FALSE)

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
