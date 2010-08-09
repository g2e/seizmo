function [x]=lind(varargin)
%LIND    Returns a linear indices matrix
%
%    Usage:    lind
%              lind(n)
%              lind(m,n) or lind([m n])
%              lind(m,n,p,...) or lind([m n p ...])
%              lind(...,classname)
%
%    Description:
%     LIND returns 1.
%
%     LIND(N) returns a NxN matrix of linear indices.
%
%     LIND(M,N) OR LIND([M N]) returns a MxN matrix of linear indices.
%
%     LIND(M,N,P,...) OR LIND([M N P ...]) returns a MxNxPx... matrix
%     of linear indices.
%
%     LIND(...,CLASSNAME) returns a linear indices matrix with class
%     CLASSNAME.
%
%    Notes:
%
%    Examples:
%     % get linear indices for some matrix 'a'
%     li=lind(size(a));
%
%    See also: ONES, ZEROS, NAN, EYE, RAND, RANDN, TRUE, FALSE

%     Version History:
%        Aug.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 12:25 GMT

% todo:

% pass to ones and then assign indices
x=ones(varargin{:});
x(:)=1:numel(x);

end
