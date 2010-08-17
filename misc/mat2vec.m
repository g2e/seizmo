function [varargout]=mat2vec(varargin)
%MAT2VEC    Converts matrices to column vectors
%
%    Usage:    vector=mat2vec(matrix)
%              [vec1,...,vecN]=mat2vec(mat1,...,matN)
%
%    Description:
%     VECTOR=MAT2VEC(MATRIX) flattens the matrix MATRIX into a column
%     vector VECTOR.  This is the functional form of VECTOR=MATRIX(:).
%
%     [VEC1,...,VECN]=MAT2VEC(MAT1,...,MATN) flattens all inputs into
%     column vectors.
%
%    Notes:
%
%    Examples:
%     % these produce the same result
%     mat=magic(5);
%     vec=mat(:)
%     vec=mat2vec(mat)
%
%    See also: COLON OPERATOR (:), NUMEL

%     Version History:
%        Aug. 14, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 14, 2010 at 13:40 GMT

% todo:

% make column vector
for i=1:nargin; varargout{i}=varargin{i}(:); end

end
