function [lgc,num,ns]=isequalnumelorscalar(varargin)
%ISEQUALNUMELORSCALAR    True if all inputs have equal numel or are scalar
%
%    Usage:    tf=isequalnumelorscalar(A,B)
%              tf=isequalnumelorscalar(A,B,C,...)
%              [tf,sz,ns]=isequalnumelorscalar(...)
%
%    Description:
%     TF=ISEQUALNUMELORSCALAR(A,B) is 1 if the two arrays have an equal
%     number of elements (that is NUMEL(A) and NUMEL(B) return the same
%     value) or are scalar.  If one is scalar the other may be any size.
%     See ISEQUALSIZEORSCALAR to require equal sized inputs.
%
%     ISEQUALNUMELORSCALAR(A,B,C,...) is 1 if all the input arguments have
%     equal number of elements or are scalar.
%
%     [TF,NUM,NS]=ISEQUALNUMELORSCALAR(...) also returns the number of
%     elements NUM and a logical array NS indicating which inputs were
%     nonscalars.  This is useful for expanding scalars (ie EXPANDSCALARS).
%
%    Notes:
%     - returns true for no input or single input cases
%
%    Examples:
%     % Will return true:
%     isequalnumelorscalar(1:10,(1:10)',19,'a')
%
%    See also: ISEQUALSIZEORSCALAR, ISEQUAL, NUMEL, EXPANDSCALARS,
%              FLATTENARRAYS

%     Version History:
%        Aug. 14, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 13:40 GMT

% todo:

% assume true
lgc=true; num=[]; ns=[];

% nobody to compare to
if(nargin==0); return; end

% single input case
if(nargin==1)
    num=numel(varargin{1});
    ns=~isscalar(varargin{1});
    return;
end

% find nonscalars
num=cellfun('prodofsize',varargin);
ns=num~=1;

% not enough nonscalars to compare
if(sum(ns)==0)
    num=1;
    return;
elseif(sum(ns)==1)
    num=numel(varargin{ns});
    return;
end

% compare nonscalars
num=unique(num(ns));
lgc=isscalar(num);

end
