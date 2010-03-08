function [varargout]=elementary_mt(n)
%ELEMENTARY_MT    Returns one of six elementary moment tensors
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

% the six elementary moment tensors
m=[ 0  0  0  1  0  0
    1 -1  0  0  0  0
    0  0  0  0  0  1
    0  0  0  0  1  0
   -1  0  1  0  0  0
    1  1  1  0  0  0];

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check n
if(~isreal(n) || any(~ismember(n,1:6)))
    error('seizmo:elementary_mt:badInput',...
        'N must be 1, 2, 3, 4, 5 or 6!');
end
n=n(:);

% output as either a single array or individually
if(nargout<=1)
    varargout{1}=m(n,:);
elseif(nargout<=6)
    m=mat2cell(m(n,:),numel(n),ones(6,1));
    [varargout{1:nargout}]=deal(m{1:nargout});
else
    error('seizmo:elementary_mt:badOutput',...
        'Incorrect number of outputs!');
end

end
