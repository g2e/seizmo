function [varargout]=seizmoverbose(varargin)
%SEIZMOVERBOSE    Turn verbose SEIZMO output on (TRUE) or off (FALSE)

% retrieve global
global SEIZMO

% get current value
if(nargout)
    try
        varargout{1}=SEIZMO.SEIZMOVERBOSE;
        if(~islogical(varargout{1}) || ~isscalar(varargout{1}))
            error('seizmo:seizmoverbose:badState',...
                'STATE of SEIZMOVERBOSE must be TRUE or FALSE!');
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
