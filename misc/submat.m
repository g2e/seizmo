function [X]=submat(X,varargin)
%SUBMAT    Returns a submatrix reduced along indicated dimensions
%
%    Usage:    Y=submat(X,DIM1,LIST1,DIM2,LIST2,...)
%
%    Description: Y=SUBMAT(X,DIM,LIST) creates a matrix Y that is
%     the matrix X reduced along dimension DIM to the indices in LIST.  If
%     DIM is a list of dimensions, LIST is used to reduce each dimension.
%
%     Y=SUBMAT(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.
%
%    Notes:
%
%    Examples:
%     Return x reduced to only the elements in index 1 of dimension 5:
%      x=submat(x,5,1)
%
%     Reduce dimensions 1 thru 3 to return a matrix that only contains
%     elements that were in row/column/page 4 or 5 for those dimensions:
%      x=submat(x,1:3,4:5)
%
%    See also: submat_eval, colon operator (:), repmat

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - fixed error message
%        Sep. 21, 2009 - updated examples, removed unnecessary brackets
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 21, 2009 at 20:45 GMT

% todo:

% CHECK VARARGIN
if(~mod(nargin,2))
    error('seizmo:submat:badNumArgs','Unpaired DIM,LIST!');
end

% DEFAULT TO ENTIRE MATRIX AND EXPAND TO MAX INPUT DIMENSION
list(1:max([ndims(X) [varargin{1:2:end}]]))={':'};

% REDUCED/REPLICATED DIMENSIONS
for i=1:2:nargin-2
    [list{varargin{i}}]=deal(varargin{i+1});
end

% SUBSET
X=X(list{:});

end
