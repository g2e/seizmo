function [li]=midx2li(sx,midx,dim)
%MIDX2LI    Translates min/max indices to linear indices
%
%    Usage:    li=midx2li(sx,midx)
%              li=midx2li(sx,midx,dim)
%
%    Description:
%     LI=MIDX2LI(SX,MIDX) takes the indices returned from either the min or
%     max function and returns linear indices.  MIDX2LI will operate on the
%     first non-singleton dimension just like min/max do without a
%     dimension argument.  SX must be the size of the matrix passed to
%     min/max (use size(x)).  MIDX is the 2nd output from min/max.
%
%     LI=MIDX2LI(SX,MIDX,DIM) specifies the dimension that min/max operated
%     along.
%
%    Notes:
%     - Special handling of min/max indices == 0 (sets linear index == 0)
%
%    Examples:
%     % Verify if midx2li can match max output for a random matrix:
%     a=rand(5,1,4,2,3);
%     [b,c]=max(a,[],1);
%     isequal(b,a(midx2li(size(a),c,1)))
%     [b,c]=max(a,[],2);
%     isequal(b,a(midx2li(size(a),c,2)))
%     [b,c]=max(a,[],3);
%     isequal(b,a(midx2li(size(a),c,3)))
%     [b,c]=max(a,[],4);
%     isequal(b,a(midx2li(size(a),c,4)))
%     [b,c]=max(a,[],5);
%     isequal(b,a(midx2li(size(a),c,5)))
%
%    See also: MAX, MIN, SORT2LI

%     Version History:
%        Feb. 15, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2011 at 16:05 GMT

% todo:

% dimension operated on
if(nargin<3 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% expand dimension list with singletons to dim
sx(end+1:dim)=1;

% linear index step size between elements along dimension
step=max([1 prod(sx(1:dim-1))]);

% linear indices of first elements along dim
li=reshape(1:prod(sx),sx);
list(1:numel(sx))={':'};
list{dim}=ones;
li=li(list{:});

% linear indices of all elements
li=li+(midx-1).*step;
li(midx==0)=0;

end
