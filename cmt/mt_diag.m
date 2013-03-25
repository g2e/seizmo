function [mt,vec]=mt_diag(varargin)
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
%     % Check diagonalization results against the eigenvalues
%     % given by the GlobalCMT catalog for 1 CMT:
%     cmt=findcmt;
%     val=mt_diag(cmt)
%     diag([cmt.eigval1 cmt.eigval2 cmt.eigval3])
%
%    See also: MT_DECOMP, MT_UNDIAG

%     Version History:
%        June  6, 2011 - initial version
%        June 11, 2011 - improved docs, now returns eigen values & vectors
%        Mar. 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check tensor format (force 3x3xN)
error(mt_check(varargin{:}));
mt=mt_change('g',varargin{:});
n=size(mt,3);

% diagonalize one by one
vec=nan(3,3,n);
for i=1:n
    % get eigenvectors & eigenvalues
    % - eigenvalue matrix is the diagonalized moment tensor
    [vec(:,:,i),mt(:,:,i)]=eig(mt(:,:,i));
end

end
