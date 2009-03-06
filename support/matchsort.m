function [y,li]=matchsort(x,i,dim)
%MATCHSORT    Replicates a sort operation using the returned permutation indices
%
%    Description:  Replicates a previous sort operation using the 
%     permutation matrix from that previous sort.  Useful for sorting
%     multiple matrices in parallel.  If no dimension argument is supplied,
%     matchsort will operate on the first non-singleton dimension just like
%     sort does.  Second output is the permutation matrix transformed into
%     the corresponding matrix of linear indices.
%
%    Usage: y=matchsort(x,i)
%           y=matchsort(x,i,dim)
%           [y,li]=matchsort(x,i,dim)
%
%    Examples:
%     
%     Sort a and then sort the elements of b,c exactly the same way.
%      b=1./a
%      c=-a
%      [d,i]=sort(a)
%      e=matchsort(b,i)
%      f=matchsort(c,i)
%      isequal(d,1./e,-f)
%
%    See also: sort

% CHECK THAT PERMUTATION MATRIX MATCHES INPUT
sx=size(x);
if(~isequal(sx,size(i)))
    error('Inputs must match in size')
end

% DIMENSION TO OPERATE ON
if(nargin<3 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% GET LINEAR INDICES FROM PERMUTATION INDICES MATRIX
li=sort2li(i,dim);

% MATCH SORT
y=x(li);

end

function [li]=sort2li(i,dim)
%SORT2LI    Transforms indices from sort to linear indices
%
%    Description:  The indices returned from the 'sort' function are
%     difficult to reimplement for sorting multiple n-d matrices in 
%     parallel.  This function transforms the indices output from sort into
%     a matrix of linear indices so that parallel sorting is simplified.
%     If no dimension argument is supplied, sort2li will operate on the 
%     first non-singleton dimension just like sort will do without a 
%     dimension argument.
%
%    Usage: li=sort2li(sort_indices)
%           li=sort2li(sort_indices,dim)
%
%    Examples:
%      
%     Sort a and then sort the elements of b exactly the same.  
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
%    Description: Y=SUBMAT_NOEVAL(X,DIM,LIST) creates a matrix Y that is
%     the matrix X reduced along dimension DIM to the indices in LIST.  If
%     DIM is a list of dimensions, LIST is used to reduce each dimension.
%
%     Y=SUBMAT_NOEVAL(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.
%
%    Usage: Y=submat_noeval(X,DIM1,LIST1,DIM2,LIST2,...)
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

