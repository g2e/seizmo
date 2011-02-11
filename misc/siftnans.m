function [x,li]=siftnans(x,dim)
%SIFTNANS    Return input with NaNs "sifted" to end of specified dimension
%
%    Usage:    y=siftnans(x)
%              y=siftnans(x,dim)
%              [y,li]=siftnans(...)
%
%    Description:
%     Y=SIFTNANS(X) sifts nan elements in X to the end of the first
%     non-singleton without sorting the non-nan elements (ie their order
%     along the dimension is preserved).
%
%     Y=SIFTNANS(X,DIM) shifts nan elements down the specified dimension.
%
%     [Y,LI]=SIFTNANS(...) returns linear indices such that Y=X(LI).
%
%    Notes:
%
%    Examples:
%     % Make a random 5x5 matrix, set some elements to nan, and sift them
%     % to the end by column then by row:
%     a=rand(5);
%     a(logical(mod(floor(a*1000),2)))=nan
%     a=siftnans(a)
%     a=siftnans(a,2)
%
%    See also: ISNAN, SORT, SIFTFUN

%     Version History:
%        Feb.  7, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% dimension to operate on
sx=size(x);
if(nargin<2 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% permute dimension to lead
x=permute(x,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);

% shift down nans
nans=isnan(x);
x(sort(~nans,1,'descend'))=x(~nans);
x(sort(nans,1))=nan;

% indices matrix if desired
if(nargout>1)
    % linear indices permuted
    li=permute(reshape(1:numel(x),sx),...
        [dim 1:dim-1 dim+1:max(numel(sx),dim)]);
    % sifted
    li([find(sort(~nans,1,'descend')); find(sort(nans,1))])=...
        li([find(~nans); find(nans)]);
    % unpermuted
    li=ipermute(li,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);
end

% permute back
x=ipermute(x,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);

end
