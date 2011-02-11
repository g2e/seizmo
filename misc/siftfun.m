function [x,li]=siftfun(x,func,dim)
%SIFTFUN    Return "sifted" input using the specified function
%
%    Usage:    y=siftfun(x,func)
%              y=siftfun(x,func,dim)
%              [y,li]=siftfun(...)
%
%    Description:
%     Y=SIFTFUN(X,FUNC) sifts elements in X identified by the function
%     underlying the function handle FUNC to the end of the first
%     non-singleton without sorting internally the sifted or residual
%     element sets.
%
%     Y=SIFTFUN(X,FUNC,DIM) shifts elements down the specified dimension.
%
%     [Y,LI]=SIFTFUN(...) returns linear indices such that Y=X(LI).
%
%    Notes:
%
%    Examples:
%     % Sift values less than .3 to bottom of an array with some nans:
%     a=rand(5);
%     a(logical(mod(floor(a*1000),2)))=nan
%     siftfun(a,@(x)x<.3)
%
%    See also: FUNCTION_HANDLE, SORT, SIFTNANS

%     Version History:
%        Feb.  7, 2011 - initial version
%        Feb. 10, 2011 - cleaned up documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% dimension to operate on
sx=size(x); nx=prod(sx);
if(nargin<3 || isempty(dim))
    dim=find(sx~=1,1);
    if(isempty(dim)); dim=1; end
end

% check func is a func
if(~isa(func,'function_handle'))
    error('misc:siftfun:badInput',...
        'FUNC must be a function handle!');
end

% permute dimension to lead
x=permute(x,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);

% sift x values
sift=func(x);

% require element by element logical output
if(~islogical(sift) || ~isequal(numel(sift),nx))
    error('misc:siftfun:badInput',...
        'FUNC must generate a logical array equal in size to X!');
end

% shift sifted
x([find(sort(~sift,1,'descend')); find(sort(sift,1))])=...
    x([find(~sift); find(sift)]);

% indices matrix if desired
if(nargout>1)
    % linear indices permuted
    li=permute(reshape(1:nx,sx),[dim 1:dim-1 dim+1:max(numel(sx),dim)]);
    % sifted
    li([find(sort(~sift,1,'descend')); find(sort(sift,1))])=...
        li([find(~sift); find(sift)]);
    % unpermuted
    li=ipermute(li,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);
end

% permute back
x=ipermute(x,[dim 1:dim-1 dim+1:max(numel(sx),dim)]);

end
