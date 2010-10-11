function [varargout]=versioninfo_cache(varargin)
%VERSIONINFO_CACHE    Check/Change VERSIONINFO caching ON=TRUE / OFF=FALSE
global SEIZMO
if(nargout)
    try
        varargout{1}=SEIZMO.VERSIONINFO.USECACHE;
    catch
        varargout{1}=false;
    end
end
if(nargin)
    SEIZMO.VERSIONINFO.USECACHE=varargin{1};
end
end
