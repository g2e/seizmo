function [varargout]=seizmocheck_state(varargin)
%SEIZMOCHECK_STATE    Check/Change if SEIZMOCHECK is ON=TRUE / OFF=FALSE
global SEIZMO
if(nargout)
    try
        varargout{1}=SEIZMO.SEIZMOCHECK.ON;
    catch
        varargout{1}=true;
    end
end
if(nargin)
    SEIZMO.SEIZMOCHECK.ON=varargin{1};
end
end
