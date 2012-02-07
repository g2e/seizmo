function [v]=mt_g2v(g)
%MT_G2V    Convert moment tensor from 3x3xN to Nx6
%
%    Usage:    mt=mt_g2v(momten)
%
%    Description:
%     MT=MT_G2V(MOMTEN) converts moment tensors in tensor form (3x3xN) to a
%     compact lower triangle form (Nx6).  This reduces memory burden (uses
%     2/3rds the memory) without losing any information assuming the tensor
%     is symmetric.
%
%    Notes:
%     - Throws an error if the tensor is not symmetric.
%
%    Examples:
%     % Compare tensor and compact forms:
%     momten=mt_v2g(elementary_mt(1:6))
%     mt_g2v(momten)
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_V2G, MT_G2C, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%        Feb. 11, 2011 - mass nargchk fix
%        June  1, 2011 - doc update
%        Feb.  7, 2012 - error if asymmetric, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check v
sz=size(g);
if(~isreal(g) || ~isequal(sz(1:2),[3 3]))
    error('seizmo:mt_g2v:badInput',...
        'Input must be a real-valued 3x3xN array!');
elseif(~isequal(g,permute(g,[2 1 3])))
    error('seizmo:mt_g2v:badInput',...
        'Cannot convert asymmetric tensors!');
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
