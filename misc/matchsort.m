function [x,li]=matchsort(x,i,dim)
%MATCHSORT    Replicates a sort operation using the returned sort indices
%
%    Usage: y=matchsort(x,sort_idx)
%           y=matchsort(x,sort_idx,dim)
%           [y,li]=matchsort(...)
%
%    Description:
%     Y=MATCHSORT(X,SORT_IDX) replicates a previous sort operation using
%     the permutation matrix from that previous sort.  This is useful for
%     sorting multiple equal-sized matrices in parallel.  Operates on the
%     first non-singleton dimension just like sort does.
%
%     Y=MATCHSORT(X,SORT_IDX,DIM) allows specifying the dimension that the
%     sort operation occured on.
%
%     [Y,LI]=MATCHSORT(...) also returns the permutation matrix transformed
%     into the corresponding matrix of linear indices for really quick
%     implementation to the sorted arrangement.  See the example below to
%     see how to use this.
%
%    Notes:
%
%    Examples:
%     % So lets say you have 5 n-dimensional same-sized matrices and you
%     % wanted to rearrange the elements in all the matrices the exact same
%     % way.  Lets say that rearrangement is based on a sorting on one of
%     % the matrices.  Sort doesn't really make this easy, so normally you
%     % would combine all the matrices into a 2D array and do a sort rows
%     % operation -- or you could use matchsort:
%     [m,i]=sort(a,4);
%     [n,li]=matchsort(b,i,4);
%     o=c(li); p=d(li); q=e(li);
%
%    See also: SORT, UNSORT, SORT2LI

%     Version History:
%        Sep.  8, 2009 - doc cleanup, dropped submat, sort2li subfunctions
%        Sep. 21, 2009 - minor doc update, changed arg y to x
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 20:35 GMT

% todo:

% CHECK THAT PERMUTATION MATRIX MATCHES INPUT
sx=size(x);
if(~isequal(sx,size(i)))
    error('seizmo:matchsort:badInput',...
        'X and SORT_IDX must match in size!');
end

% DIMENSION TO OPERATE ON
if(nargin<3 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% GET LINEAR INDICES FROM PERMUTATION INDICES MATRIX
li=sort2li(i,dim);

% MATCH SORT
x=x(li);

end
