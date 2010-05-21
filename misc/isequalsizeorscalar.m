function [lgc,sz,ns]=isequalsizeorscalar(varargin)
%ISEQUALSIZEORSCALAR    True if all input arrays are equal size or scalar
%
%    Usage:    isequalsizeorscalar(A,B)
%              isequalsizeorscalar(A,B,C,...)
%              [tf,sz,ns]=isequalsizeorscalar(...)
%
%    Description:  ISEQUALSIZEORSCALAR(A,B) is 1 if the two arrays are the
%     same size or are scalar.  If one is scalar the other may be any size.
%
%     ISEQUALSIZEORSCALAR(A,B,C,...) is 1 if all the input arguments have
%     equal size or are scalar.
%
%     [TF,SZ,NS]=ISEQUALSIZEORSCALAR(...) also returns the nonscalar array
%     size SZ and a logical array NS indicating which inputs were
%     nonscalars.  This is useful for expanding scalars (ie EXPANDSCALARS).
%
%    Notes:
%     - returns true for no input or single input cases
%
%    Examples:
%     Make some arrays to test:
%      A=cell(3,4,10);
%      B=nan(3,4,10);
%      C=4;
%      D=magic(5);
%
%     Will return true:
%      isequalsizeorscalar(A,B,C)
%
%     Will return false:
%      isequalsizeorscalar(A,B,C,D)
%
%    See also: ISEQUAL, SIZE, EXPANDSCALARS

%     Version History:
%        May   1, 2010 - initial version
%        May  16, 2010 - more outputs
%        May  20, 2010 - no split based on matlab version (sloooow)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  20, 2010 at 13:40 GMT

% todo:

% assume true
lgc=true; sz=[]; ns=[];

% nobody to compare to
if(nargin==0); return; end

% single input case
if(nargin==1); sz=size(varargin{1}); ns=~isscalar(varargin{1}); return; end

% find nonscalars
ns=cellfun('prodofsize',varargin)~=1;

% not enough nonscalars to compare
if(sum(ns)==0)
    sz=[1 1];
    return;
elseif(sum(ns)==1)
    sz=size(varargin{ns});
    return;
end

% compare nonscalars
% - this breaks with early cellfun
%if(verLessThan('matlab', '7.1'))
%    nsi=find(ns);
%    sz=size(varargin{nsi(1)});
%    for i=nsi(2:end)
%        if(~isequal(sz,size(varargin{i}))); lgc=false; return; end
%    end
%else % 7.1+
    sizes=cellfun(@size,varargin(ns),'UniformOutput',false);
    lgc=isequal(sizes{:});
    sz=sizes{1};
%end

end
