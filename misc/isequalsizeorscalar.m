function [lgc]=isequalsizeorscalar(varargin)
%ISEQUALSIZEORSCALAR    True if all input arrays are equal size or scalar
%
%    Usage:    isequalsizeorscalar(A,B)
%              isequalsizeorscalar(A,B,C,...)
%
%    Description:  ISEQUALSIZEORSCALAR(A,B) is 1 if the two arrays are the
%     same size or are scalar.  If one is scalar the other may be any size.
%
%     ISEQUALSIZEORSCALAR(A,B,C,...) is 1 if all the input arguments have
%     equal size or are scalar.
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
%    See also: ISEQUAL, SIZE

%     Version History:
%        May   1, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   1, 2010 at 09:25 GMT

% todo:

% assume true
lgc=true;

% nobody to compare to
if(nargin<2); return; end

% find nonscalars
ns=cellfun('prodofsize',varargin)~=1;

% not enough nonscalars to compare
if(sum(ns)<2); return; end

% compare nonscalars
% - this breaks with early cellfun
if(verLessThan('matlab', '7.1'))
    nsi=find(ns);
    sz=size(varargin{nsi(1)});
    for i=nsi(2:end)
        if(~isequal(sz,size(varargin{i}))); lgc=false; return; end
    end
else % 7.1+
    sizes=cellfun(@size,varargin(ns),'UniformOutput',false);
    lgc=isequal(sizes{:});
end

end
