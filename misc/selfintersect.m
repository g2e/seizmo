function [p,sidx]=selfintersect(poly,y)
%SELFINTERSECT    Self intersection points & segments of a 2D polygon
%
%    Usage:    [p,segidx]=selfintersect(xy)
%              [p,segidx]=selfintersect(x,y)
%
%    Description:
%     [P,SEGIDX]=SELFINTERSECT(XY) locates self-intersection points of a 2D
%     polygon given by XY.  XY is a Nx2 or 2xN array of x & y positions for
%     the N vertices of the polygon.  P is a 2xM array of the intersection
%     points found.  SEGIDX is a 1xM cell-array of polygon segment indices
%     for each intersection point (each cell has 2 indices).
%
%     [P,SEGIDX]=SELFINTERSECT(X,Y) allows specifying the x & y positions
%     of the polygon vertices independently.  Scalar expansion is not
%     supported (for obvious reasons).
%
%    Notes:
%     - Segment indices can be converted to indices fairly easily:
%        segment   1 = vectices   1 &  2
%        segment   2 = vertices   2 &  3
%        ...
%        segment N-1 = vertices N-1 &  N
%        segment   N = vertices   N &  1
%
%    Examples:
%     % Find and plot intersection info of a 10pt polygon:
%     x=rand(1,10); x=[x x(1)];
%     y=rand(1,10); y=[y y(1)];
%     [p,sidx]=selfintersect(x,y);
%     fh=figure; ax=axes('parent',fh);
%     plot(ax,x,y);
%     hold(ax,'on');
%     plot(ax,p(1,:),p(2,:),'r*');
%     plot(ax,x([cell2mat(sidx); cell2mat(sidx)+1]),...
%             y([cell2mat(sidx); cell2mat(sidx)+1]),'r');
%     hold(ax,'off');
%
%    See also: POLYAREA, INPOLYGON

%     Version History:
%        Jan. 23, 2011 - initial version
%        Feb. 10, 2011 - revised to return segment indices, major doc
%                        update, code refactoring (drop buffer code)
%
%     Written by Bruno Luong <brunoluong@????.???>
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));
if(nargin==2); poly=[poly(:) y(:)].';
else if(size(poly,2)==2 && size(poly,2)~=2); poly=poly.'; end;
end

% wrap around if needed
if(~isequal(poly(:,end),poly(:,1))); poly(:,end+1)=poly(:,1); end

% Loop over segment of poly
p=zeros(2,0); sidx={};
for n=2:size(poly,2)-2
    seg=poly(:,n-1:n);
    subpoly=poly(:,n+1:end);
    if(n==2); subpoly(:,end)=[]; end
    [xp,idx]=segxpoly(seg,subpoly);
    nidx=numel(idx);
    if(~nidx); continue; end
    idx=mat2cell(...
        reshape([(n-1)*ones(1,nidx); n+idx],1,2*nidx),1,2*ones(1,nidx));
    p=[p xp]; %#ok<AGROW>
    sidx=[sidx idx]; %#ok<AGROW>
end

end


function [p,idx]=segxpoly(seg,poly)
% function [p,idx]=segxpoly(seg,poly)
% Check if a line segment seg intersects with a polygon poly.
% INPUTS:
% seg is (2 x 2) where
% seg(:,1) is the first point
% seg(:,2) is the the second point of the segment.
% poly is (2 x n) array, each column is a vertices
% OUTPUT
% p is (2 x m) array, each column is an intersecting point
% idx is (1 x m) vector of intersected polygon segment indices (start pt)
%
% Author: Bruno Luong <brunoluong@????.???>
% History:
% Original 20-May-2010
% Edited 10-Feb-2011 (gge) return indices too

% Translate so that first point is origin
a=seg(:,1);
M=bsxfun(@minus,poly,a);
b=seg(:,2)-a;
% Check if the points are on the left/right side
x=[b(2) -b(1)]*M;
sx=sign(x);
% x-coordinates has opposite signs
idx=sx(1:end-1).*sx(2:end)<=0;
if(any(idx))
    idx=find(idx);
    % cross point to the y-axis (along the segment)
    x1=x(idx);
    x2=x(idx+1);
    d=b.'/(b(1)^2+b(2)^2);
    y1=d*M(:,idx);
    y2=d*M(:,idx+1);
    dx=x2-x1;
    % We won't bother with the degenerate case of dx=0 and x1=0
    y=(y1.*x2-y2.*x1)./dx;
    % Check if the cross point is inside the segment
    idx2=y>=0 & y<1;
    if(any(idx2))
        p=bsxfun(@plus,a,b*y(idx2));
        idx=idx(idx2);
    else
        p=zeros(2,0);
        idx=zeros(1,0);
    end
else
    p=zeros(2,0);
    idx=zeros(1,0);
end

end
