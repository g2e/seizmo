function [uti]=sub2uti(nrows,i,j)
%SUB2UTI    Square matrix upper triangle linear indices from subscripts
%
%    Usage:    idx=sub2uti(i,j)
%
%    Description:
%     IDX=SUB2UTI(I,J) converts row/column subscripts I & J to upper
%     triangle indices IDX.  The number of rows in the corresponding matrix
%     does not need to be given.  Upper triangle indices increase similar
%     to linear indices, but the lower triangle and the diagonal are
%     ignored so that for a 5x5 matrix the indices proceed as follows:
%
%      \j  1  2  3  4  5
%     i \
%     1    -  1  2  4  7
%     2    -  -  3  5  8
%     3    -  -  -  6  9
%     4    -  -  -  - 10
%     5    -  -  -  -  -
%
%    Notes:
%
%    Examples:
%     % Say you have a dissimilarity vector (in this case, it corresponds
%     % to the upper triangle of a 400x400 dissimilarity matrix) and you
%     % wanted to know the dissimilarity between thingy 74 and 233:
%     idx=sub2uti(400,74,233)
%     dissim=my_dissim_vector(idx)
%
%    See also: SUB2LTI, UTI2SUB, SUB2LTI

%     Version History:
%        Sep.  8, 2009 - added documentation
%        Oct. 13, 2009 - added checks, updated documentation
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 21:40 GMT

% todo:

% checks
if(isempty(nrows) || ~isnumeric(nrows) ...
        || any(nrows~=fix(nrows)) || numel(nrows)>2)
    error('seizmo:sub2uti:badInput',...
        'NROWS must be a scalar number or 1x2 array of [NROWS NCOLS]!');
elseif(~isnumeric(i) || ~isnumeric(j) ...
        || any(i~=fix(i)) || any(j~=fix(j)) ...
        || (~isscalar(i) && ~isscalar(j) && numel(i)~=numel(j)))
    error('seizmo:sub2uti:badInput',...
        ['I and J must be scalar integers or be\n' ...
        'arrays with the same number of integers!']);
end

% get first dimension if more than 1
nrows=nrows(1);

% subscripts should be in upper triangle
if(any(i>=j))
    error('seizmo:sub2uti:badInput','Indices outside upper triangle!');
end

% get upper triangle indices
j=j(:)-1;
k=[0 1:max(j)];
k(k>nrows)=nrows;
k=cumsum(k);
uti=k(j).'+i(:);

end
