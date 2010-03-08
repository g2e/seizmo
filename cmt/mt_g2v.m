function [v]=mt_g2v(g)
%MT_G2V    Convert moment tensor from 3x3xN to Nx6
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
sz=size(g);
if(~isreal(g) || ~isequal(sz(1:2),[3 3]))
    error('seizmo:mt_g2v:badInput',...
        'G must be a real-valued 3x3xN array!');
end

% convert
nd=numel(sz);
if(nd==2); sz(3)=1; end
v=nan([6 sz(3:end)]);
v([1 4 5],:)=g(1:3,1,:);
v([2 6],:)=g(2:3,2,:);
v(3,:)=g(3,3,:);
v=permute(v,[2 1 3:nd-1]);

end
