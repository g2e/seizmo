function [varargout]=mt2sdr(varargin)
%MT2SDR    Convert moment tensor to strike-dip-rake
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

% what has to be done
% - decompose to a pure double couple
%   - how?
%     - diagonalize (eig)
%     - remove any isotropic component
%     - break into major and minor double couple
%     - use major double couple
% - now go from pure double couple to focal planes
%   - how?
%     - eigenvectors give tpb
%       - t is the eigenvector for the positive eigenvalue
%       - p is the eigenvector for the negative eigenvalue
%     - convert tpb to nd
%       - n=t+p
%       - d=t-p
%     - convert nd to sdr
%       - ...

end
