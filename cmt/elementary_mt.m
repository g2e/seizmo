function [varargout]=elementary_mt(n)
%ELEMENTARY_MT    Returns one of six elementary moment tensors
%
%    Usage:    mt=elementary_mt(idx)
%
%    Description: MT=ELEMENTARY_MT(IDX) returns elementary moment tensors
%     used in moment tensor inversion.  There are 6 different elementary
%     moment tensors: 1 & 2 are strike-slip faults, 3 & 4 are dip-slip
%     faults on vertical planes, 5 is a 45deg dip-slip fault, & 6 is an
%     isotropic source (explosion).  IDX should be a scalar or vector of
%     integers from 1 to 6.  MT is a Nx6 array of moment tensor components
%     where N is the number of integers in IDX.  MT is in Aki & Richard's
%     form (North, East, Up axis) giving [Mxx Myy Mzz Mxy Mxz Myz].
%
%    Notes:
%
%    Examples:
%     Show the 6 elementary focal mechanisms:
%      bb(ar2hrv(elementary_mt(1:6)),1:6,ones(1,6),0.5,0,'b');
%      axis equal off
%
%    See also: AR2HRV

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 22, 2010 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 22, 2010 at 15:05 GMT

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
