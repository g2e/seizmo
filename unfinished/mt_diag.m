function [d,v]=mt_diag(varargin)
%MT_DIAG    Returns diagonalized moment tensors
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

% todo:

if(nargin==1)
    if(~isreal(varargin{1}))
        error('seizmo:mt_diag:badInput',...
            'MOMTEN must be real-valued!');
    end
    sz=size(varargin{1});
    if(isequal([3 3],sz(1:2)))
        if(numel(sz)==2); sz(3)=1; end
        nmt=prod(sz(3:end));
        varargout{1}=varargin{1};
        varargout{2}=nan(nmt,1);
        for i=1:nmt
            varargout{2}(i)=sqrt(...
                trace(varargout{1}(:,:,i)*varargout{1}(:,:,i))/2);
            varargout{1}(:,:,i)=varargout{1}(:,:,i)/varargout{2}(i);
        end
    elseif(isequal(6,sz(2)))
        nmt=prod(sz([1 3:end]));
        varargout{1}=mt_v2g(varargin{1});
        varargout{2}=nan(nmt,1);
        for i=1:nmt
            varargout{2}(i)=sqrt(...
                trace(varargout{1}(:,:,i)*varargout{1}(:,:,i))/2);
            varargout{1}(:,:,i)=varargout{1}(:,:,i)/varargout{2}(i);
        end
        varargout{1}=mt_g2v(varargout{1});
    else
        error('seizmo:mt_diag:badInput',...
            'MOMTEN inproper size!');
    end
elseif(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:mt_norm:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    varargout{1}=mt_62v(varargin{:});
    sz=size(varargout{1});
    nmt=prod(sz([1 3:end]));
    varargout{1}=mt_v2g(varargout{1});
    varargout{7}=nan(nmt,1);
    for i=1:nmt
        varargout{7}(i)=sqrt(...
            trace(varargout{1}(:,:,i)*varargout{1}(:,:,i))/2);
        varargout{1}(:,:,i)=varargout{1}(:,:,i)/varargout{7}(i);
    end
    varargout{1}=mt_g2v(varargout{1});
    [varargout{1:6}]=mt_v26(varargout{1});
else
    error('seizmo:mt_norm:badNumInput',...
        'Incorrect number of inputs!');
end

end