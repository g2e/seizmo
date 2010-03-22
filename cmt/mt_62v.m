function [varargout]=mt_62v(varargin)
%MT_62V    Converts moment tensor from 6 Nx1 vectors to Nx6 array
%
%    Usage:    mt=mt_62v(M11,M22,M33,M12,M13,M23)
%
%    Description: MT=MT_62V(M11,M22,M33,M12,M13,M23) converts a moment
%     tensor given as separate component arrays to a combined format that
%     is just the 6 components concatenated along dimension 2.  This gives
%     an Nx6 array that is the most compact format.
%
%    Notes:
%
%    Examples:
%     Combine Harvard moment tensor components to a combined array:
%      mt=mt_62v(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp);
%
%    See also: MT_V26, MT_V2G, MT_G2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 21, 2010 at 11:25 GMT

% todo:

if(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:mt_62v:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    % expand scalars
    n=cellfun('prodofsize',varargin);
    sz=size(varargin{find(n==max(n),1,'first')});
    for i=1:6
        if(n(i)==1)
            varargin{i}=varargin{i}(ones(sz));
        else
            if(~isequal(size(varargin{i}),sz))
                error('seizmo:mt_62v:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    varargout{1}=cat(2,varargin{:});
else
    error('seizmo:mt_62v:badNumInput',...
        'Incorrect number of inputs!');
end

end
