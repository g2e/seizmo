function [g]=mt_v2g(v)
%MT_V2G    Convert moment tensor from Nx6 to 3x3xN
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

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check v
sz=size(v);
if(~isreal(v) || sz(2)~=6)
    error('seizmo:mt_v2g:badInput',...
        'V must be a real-valued Nx6 array!');
end

% convert
g=nan([3 3 sz([1 3:end])]);
v=permute(v,[2 1 3:numel(sz)]);
v=v([1 4 5 4 2 6 5 6 3],:,:);
g(:)=v(:);

end
