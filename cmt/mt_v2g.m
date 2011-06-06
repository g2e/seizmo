function [g]=mt_v2g(v)
%MT_V2G    Convert moment tensor from Nx6 to 3x3xN
%
%    Usage:    momten=mt_v2g(mt)
%
%    Description:
%     MOMTEN=MT_V2G(MT) converts moment tensors from the compact Nx6 form
%     to 3x3xN tensor form.  This requires more memory (50% more) but is
%     the appropriate mathematical form.  It is useful for operations that
%     operate on the actual tensor (diagonalization etc).
%
%    Notes:
%
%    Examples:
%     % Compare tensor and compact forms:
%     momten=mt_v2g(elementary_mt(1:6))
%     mt_g2v(momten)
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_G2V, MT_G2C, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%        Feb. 11, 2011 - mass nargchk fix
%        June  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check v
sz=size(v);
if(~isreal(v) || sz(2)~=6)
    error('seizmo:mt_v2g:badInput',...
        'Input must be a real-valued Nx6 array!');
end

% convert
g=nan([3 3 sz([1 3:end])]);
v=permute(v,[2 1 3:numel(sz)]);
v=v([1 4 5 4 2 6 5 6 3],:,:);
g(:)=v(:);

end
