function [x,li]=unsort(x,i,dim)
%UNSORT    Undoes a sort operation using the returned sort indices
%
%    Usage:    y=unsort(x,sort_idx)
%              y=unsort(x,sort_idx,dim)
%              [y,li]=unsort(...)
%
%    Description:
%     Y=UNSORT(X,SORT_IDX) undoes a previous sort operation using the
%     permutation matrix from that previous sort.  Really this just calls
%     SORT2LI to get the sort indices into linear indices so unsorting is
%     as simple as X(LI)=X, where LI are the linear indices.  See the third
%     usage form and the example for more on linear indices.
%
%     Y=UNSORT(X,SORT_IDX,DIM) specifies the dimension that the previous
%     sort operation worked on.  If you specified DIM with sort, then DIM
%     should be specified for unsort too!
%
%     [Y,LI]=UNSORT(...) returns the linear indices, which are much more
%     convenient for re-sorting and unsorting.  See the example for usage.
%
%    Notes:
%
%    Examples:
%     % Sort and unsort a matrix:
%     [x,i]=sort(x)
%     [x,li]=unsort(x,i)
%     
%     % Now re-sort and re-unsort using the linear indices:
%     x=x(li); % sorted
%     x(li)=x; % unsorted
%
%    See also: SORT, MATCHSORT, SORT2LI

%     Version History:
%        Sep. 21, 2009 - initial version
%        Feb. 16, 2010 - fixed usage section in docs
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 01:25 GMT

% todo:

% CHECK THAT PERMUTATION MATRIX MATCHES INPUT
sx=size(x);
if(~isequal(sx,size(i)))
    error('seizmo:unsort:badInput',...
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
x(li)=x;

end
