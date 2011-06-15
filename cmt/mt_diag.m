function [mt,vec]=mt_diag(mt)
%MT_DIAG    Returns diagonalized moment tensors & principal axes
%
%    Usage:    [val,vec]=mt_diag(mt)
%
%    Description:
%     [VAL,VEC]=MT_DIAG(MT) diagonalizes the moment tensor(s) in MT using
%     an eigendecomposition.  The eigenvalues are returned down the
%     diagonal of VAL (a 3x3xN matrix) and the eigenvectors are in VEC
%     (also a 3x3xN matrix).  The columns going from 1 to 3 give the P, B,
%     & T principal axes in that order.  MT must be a 3x3xN or Nx6 moment
%     tensor array where N is the number of moment tensors in MT.  This
%     operation is useful for investigating the nature of the moment tensor
%     source (double couple, CLVD, & isotropic components).
%
%    Notes:
%     - Principal Axes:
%       T-axis - largest eigenvalue - center of compressional quadrant
%       P-axis - smallest eigenvalue - center of dilatational quadrant
%       B-axis - middle eigenvalue - fault & auxiliary plane intersection
%                (also known as the Null or Neutral axis)
%
%    Examples:
%     % Check diagonalization results against the eigenvalues given
%     % by the GlobalCMT catalog for 1 CMT:
%     cmt=findcmt
%     val=mt_diag(mt_s2g(cmt))
%
%    See also: MT_DECOMP, MT_UNDIAG

%     Version History:
%        June  6, 2011 - initial version
%        June 11, 2011 - improved docs, now returns eigen values & vectors
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2011 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check tensor format (force 3x3xN)
mtsz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt_diag:badInput',...
        'MT must be a real-valued numeric array!');
elseif(isequal(mtsz(1:2),[3 3]) && any(numel(mtsz)==[2 3]))
    if(numel(mtsz)>2); n=mtsz(3); else n=1; end
elseif(mtsz(2)==6 && numel(mtsz)==2)
    mt=mt_v2g(mt); % convert from Nx6 to 3x3xN
    n=mtsz(1);
else
    error('seizmo:mt_diag:badInput',...
        'MT must be a harvard moment tensor array as 3x3xN or Nx6!');
end

% diagonalize one by one
vec=nan(3,3,n);
for i=1:n
    % get eigenvectors & eigenvalues
    % - eigenvalue matrix is the diagonalized moment tensor
    [vec(:,:,i),mt(:,:,i)]=eig(mt(:,:,i));
end

end
