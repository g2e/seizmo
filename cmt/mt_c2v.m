function [varargout]=mt_c2v(varargin)
%MT_C2V    Converts from 6 Nx1 component vectors to Nx6 moment tensor array
%
%    Usage:    mt=mt_c2v(M11,M22,M33,M12,M13,M23)
%
%    Description:
%     MT=MT_C2V(M11,M22,M33,M12,M13,M23) converts a moment tensor given as
%     separate component column vectors to a combined format that is just
%     the 6 components concatenated as separate columns in a Nx6 array.
%
%    Notes:
%
%    Examples:
%     % Combine Harvard moment tensor components to a combined array:
%     mt=mt_62v(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp);
%
%    See also: MT_C2G, MT_V2C, MT_V2G, MT_G2V, MT_G2C, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%        June  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2011 at 11:25 GMT

% todo:

if(nargin==6)
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:mt_c2v:badInput',...
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
                error('seizmo:mt_c2v:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    varargout{1}=cat(2,varargin{:});
else
    error('seizmo:mt_c2v:badNumInput',...
        'Incorrect number of inputs!');
end

end
