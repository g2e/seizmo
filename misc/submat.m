function [X]=submat(X,varargin)
%SUBMAT    Returns a submatrix reduced along indicated dimensions
%
%    Usage:    Y=submat(X,DIM1,LIST1,DIM2,LIST2,...)
%
%    Description:
%     Y=SUBMAT(X,DIM,LIST) creates a matrix Y that is the matrix X reduced
%     along dimension DIM to the indices in LIST.  If DIM is a list of
%     dimensions, LIST is used to reduce each dimension.  DIM may be ':' to
%     indicate all dimensions of X (up to the maximum non-singleton).
%
%     Y=SUBMAT(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.  Later inputs take indexing
%     preference if a dimension is indexed more than once.
%
%    Notes:
%
%    Examples:
%     % Return x reduced to only the elements in index 1 of dimension 5:
%     x=submat(x,5,1)
%
%     % Reduce dimensions 1 thru 3 to return a matrix that only contains
%     % elements that were in row/column/page 4 or 5 for those dimensions:
%     x=submat(x,1:3,4:5)
%
%     % Reduce to elements in the 3rd row, 4th column:
%     x=submat(x,1,3,2,4)
%
%     % Return the element with all subscripts of 2:
%     x=submat(x,':',2)
%
%    See also: SUBMAT_EVAL, COLON OPERATOR (:), REPMAT

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - fixed error message
%        Sep. 21, 2009 - updated examples, removed unnecessary brackets
%        Aug.  9, 2010 - ':' dimension now functions, 15% slower though
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 20:45 GMT

% todo:

% CHECK VARARGIN
if(~mod(nargin,2))
    error('misc:submat:badNumArgs','Unpaired DIMENSION,INDICES!');
end

% FIX ':' DIMENSION
varargin(2.*find(strcmp(':',varargin(1:2:end)))-1)={1:ndims(X)};

% DEFAULT TO ENTIRE MATRIX AND EXPAND TO MAX INPUT DIMENSION
list(1:max([ndims(X) [varargin{1:2:end}]]))={':'};

% REDUCED/REPLICATED DIMENSIONS
for i=1:2:nargin-2
    [list{varargin{i}}]=deal(varargin{i+1});
end

% SUBSET
X=X(list{:});

end
