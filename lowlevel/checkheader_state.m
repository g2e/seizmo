function [varargout]=checkheader_state(varargin)
%CHECKHEADER_STATE    Check/Change if CHECKHEADER is ON=TRUE / OFF=FALSE
global SEIZMO
if(nargout)
    try
        varargout{1}=SEIZMO.CHECKHEADER.ON;
    catch
        varargout{1}=true;
    end
end
if(nargin)
    SEIZMO.CHECKHEADER.ON=varargin{1};
end
end
