function [g]=ndsquareform(v,flag)
%NDSQUAREFORM    Reshapes a multi-page distance matrix between square and triangle vector forms
%
% Usage: g=ndsquareform(v,flag)
%
%        Optional flag = TRUE (DEFAULT)
%                               -> gives row vector for grid-->vector
%                               -> gives symmetric grid for vector-->grid
%                 flag = FALSE 
%                               -> gives column vector for grid-->vector
%                               -> gives anti-symmetric grid for vector-->grid
%
% Note: 1x1x... is always assumed a vector, so 1x1x... squares are NOT
%       supported.
%
% Example:
%      v=repmat(1:10,[1 1 2 2])
%      g=ndsquareform(v)
%      g(:,:,1,1)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      g(:,:,2,1)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      g(:,:,1,2)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      g(:,:,2,2)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      v=squareform(g)
%
% See also: squareform, tril, triu, permute

% flag
if(nargin<2 || isempty(flag)); flag=true; end

% get size
sz=size(v);
ndim=length(sz);
pages=max([1 prod(sz(3:end))]);

% vector 2 grid
if(any(sz(1:2)==1))
    len=prod(sz(1:2));
    n=ceil(sqrt(2*len));
    if(len~=(n^2-n)/2); error('vector size bad'); end
    g=zeros([n n sz(3:end)]);
    lt=tril(true(n),-1); % lower triangle to preserve indices
    g(lt(:,:,ones(1,pages)))=v(:);
    if(flag); g=g+permute(g,[2 1 3:ndim]);
    else g=permute(g,[2 1 3:ndim])-g; end
% grid 2 vector
else
    if(sz(1)~=sz(2)); error('grid not square'); end
    if(~isequal(v,permute(v,[2 1 3:ndim])) ...
            && ~isequal(v,-permute(v,[2 1 3:ndim])))
        error('grid not (anti-)symmetric'); 
    end
    v=permute(v,[2 1 3:ndim]); % swap to preserve indices
    lt=tril(true(sz(1)),-1);
    g=v(lt(:,:,ones(1,pages)));
    if(flag); g=permute(g,[2 1 3:ndim]); end
end

end
