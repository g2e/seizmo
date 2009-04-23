function [li]=sort2li(i,dim)
%SORT2LI    Transforms permutation indices from sort to linear indices
%
%    Usage: li=sort2li(sort_indices)
%           li=sort2li(sort_indices,dim)
%
%    Description:  The indices returned from the 'sort' function are
%     difficult to reimplement for sorting multiple n-d matrices in 
%     parallel.  This function transforms the indices output from sort into
%     a matrix of linear indices so that parallel sorting is simplified.
%     If no dimension argument is supplied, sort2li will operate on the 
%     first non-singleton dimension just like sort will do without a 
%     dimension argument.
%
%    Examples:
%      
%     Sort a and then sort the elements of b exactly the same way.  
%      b=1./a;
%      [c,i]=sort(a)
%      d=b(sort2li(i))
%      isequal(c,1./d)
%
%    See also: sort

% DIMENSION TO OPERATE ON
sx=size(i);
if(nargin<2 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% EXPAND DIMENSION LIST WITH SINGLETONS TO DIM
sx(end+1:dim)=1;

% LINEAR INDEX STEP SIZE BETWEEN ELEMENTS ALONG DIMENSION
step=max([1 prod(sx(1:dim-1))]);

% LINEAR INDICES OF FIRST ELEMENTS ALONG DIM
li=reshape(1:numel(i),sx);
li=submat_noeval(li,dim,ones(1,sx(dim)));

% LINEAR INDICES OF ALL ELEMENTS
li=reshape(li+(i-1).*step,sx);

end

function [X]=submat_noeval(X,varargin)
%SUBMAT_NOEVAL    Returns a submatrix reduced along indicated dimensions
%
%    Usage: Y=submat_noeval(X,DIM1,LIST1,DIM2,LIST2,...)
%
%    Description: Y=SUBMAT_NOEVAL(X,DIM,LIST) creates a matrix Y that is
%     the matrix X reduced along dimension DIM to the indices in LIST.  If
%     DIM is a list of dimensions, LIST is used to reduce each dimension.
%
%     Y=SUBMAT_NOEVAL(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.
%
%    Examples:
%      Return x reduced to only the elements in index 1 of dimension 5:
%      x=submat_noeval(x,5,1)
%
%    See also: colon operator (:), repmat

% CHECK VARARGIN
if(~mod(nargin,2))
    error('dimension argument must be followed by indices argument');
end

% DEFAULT TO ENTIRE MATRIX AND EXPAND TO MAX INPUT DIMENSION
list(1:max([ndims(X) [varargin{1:2:end}]]))={':'};

% REDUCED/REPLICATED DIMENSIONS
for i=1:2:nargin-2
    [list{[varargin{i}]}]=deal(varargin{i+1});
end

% SUBSET
X=X(list{:});

end
