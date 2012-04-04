function [li]=sort2li(i,dim)
%SORT2LI    Transforms permutation indices from sort to linear indices
%
%    Usage:    li=sort2li(sort_idx)
%              li=sort2li(sort_idx,dim)
%
%    Description:
%     LI=SORT2LI(SORT_IDX) takes the indices returned from the 'sort'
%     function, which are difficult to reimplement for sorting multiple
%     n-dimensional matrices in parallel, and returns linear indices.
%     SORT2LI will operate on the first non-singleton dimension just like
%     sort does without a dimension argument.
%
%     LI=SORT2LI(SORT_IDX,DIM) allows specifying the dimension that sort
%     operated on.
%
%    Notes:
%
%    Examples:
%     % Sort a and then sort the elements of b exactly the same way.  
%     b=1./a;
%     [c,i]=sort(a)
%     d=b(sort2li(i))
%     isequal(c,1./d)
%
%    See also: SORT, SUBMAT, MATCHSORT, UNSORT, MIDX2LI

%     Version History:
%        Sep.  8, 2009 - doc cleanup, dropped submat subfunction
%        Sep. 21, 2009 - dropped submat call (now uses only built-ins)
%        Feb. 15, 2011 - added midx2li to See also section
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 21:00 GMT

% todo:

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
list(1:numel(sx))={':'};
list{dim}=ones(1,sx(dim));
li=li(list{:});

% LINEAR INDICES OF ALL ELEMENTS
li=reshape(li+(i-1).*step,sx);

end
