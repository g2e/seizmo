function [g]=mt_v2g(v)
%MT_V2G    Convert moment tensor from Nx6 to 3x3xN
%
%    Usage:    momten=mt_v2g(mt)
%
%    Description: MOMTEN=MT_V2G(MT) converts moment tensors from the
%     compact Nx6 form to 3x3xN tensor form.  This requires more memory
%     (50% more) but is the standard.  It is useful for operations that
%     operate on the actual tensor (diagonalization etc).
%
%    Notes:
%
%    Examples:
%     Compare tensor and compact forms:
%      momten=mt_v2g(elementary_mt(1:6));
%      momten(:,:,:)
%      mt_g2v(momten)
%
%    See also: MT_62V, MT_V26, MT_G2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 21, 2010 at 11:45 GMT

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
