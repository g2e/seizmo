function [varargout]=mt_c2g(varargin)
%MT_C2G    Converts from 6 Nx1 component vectors to a 3x3xN moment tensor
%
%    Usage:    momten=mt_c2g(M11,M22,M33,M12,M13,M23)
%
%    Description:
%     MOMTEN=MT_C2G(M11,M22,M33,M12,M13,M23) converts moment tensors given
%     as separate component column vectors into the moment tensor format
%     where each tensor is a separate "page" (3rd dimension) in a 3x3xN
%     array.  Components must be scalar or equal sized column vectors.
%
%    Notes:
%
%    Examples:
%     % Combine Harvard moment tensor components to a combined array:
%     momten=mt_62g(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp);
%
%    See also: MT_C2V, MT_V2C, MT_V2G, MT_G2V, MT_G2C, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        June  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2011 at 11:25 GMT

% todo:

if(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:mt_c2g:badInput',...
            'All inputs must be real-valued Nx1 column vectors!');
    end
    % expand scalars
    n=cellfun('prodofsize',varargin);
    sz=size(varargin{find(n==max(n),1,'first')});
    for i=1:6
        if(n(i)==1)
            varargin{i}=varargin{i}(ones(sz));
        else
            if(~isequal(size(varargin{i}),sz))
                error('seizmo:mt_c2g:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    varargout{1}=cat(3,...
        [varargin{[1 4 5]}],[varargin{[4 2 6]}],[varargin{[5 6 3]}]);
    varargout{1}=permute(varargout{1},[2 3 1]);
else
    error('seizmo:mt_c2g:badNumInput',...
        'Incorrect number of inputs!');
end

end
