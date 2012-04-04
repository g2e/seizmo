function [lti]=sub2lti(nrows,i,j)
%SUB2LTI    Square matrix lower triangle linear indices from subscripts
%
%    Usage:    idx=sub2lti(nrows,i,j)
%
%    Description:
%     IDX=SUB2LTI(NROWS,I,J) converts row/column subscripts I & J to lower
%     triangle indices IDX.  The number of rows in the corresponding matrix
%     needs to be given as NROWS.  Lower triangle indices increase similar
%     to linear indices, but the upper triangle and the diagonal are
%     ignored so that for a 5x5 matrix the indices proceed as follows:
%
%      \j  1  2  3  4  5
%     i \
%     1    -  -  -  -  -
%     2    1  -  -  -  -
%     3    2  5  -  -  -
%     4    3  6  8  -  -
%     5    4  7  9 10  -
%
%    Notes:
%
%    Examples:
%     % Say you have a dissimilarity vector (in this case it corresponds
%     % to the lower triangle of a 400x400 dissimilarity matrix) and you
%     % wanted to know the dissimilarity between thingy 74 and 233:
%     idx=sub2lti(400,74,233)
%     dissim=my_dissim_vector(idx)
%
%    See also: SUB2UTI, LTI2SUB, UTI2SUB

%     Version History:
%        Sep.  8, 2009 - added documentation
%        Oct. 13, 2009 - added checks, updated documentation
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 21:30 GMT

% todo:

% checks
if(isempty(nrows) || ~isnumeric(nrows) ...
        || any(nrows~=fix(nrows)) || numel(nrows)>2)
    error('seizmo:sub2lti:badInput',...
        'NROWS must be a scalar number or 1x2 array of [NROWS NCOLS]!');
elseif(~isnumeric(i) || ~isnumeric(j) ...
        || any(i~=fix(i)) || any(j~=fix(j)) ...
        || (~isscalar(i) && ~isscalar(j) && numel(i)~=numel(j)))
    error('seizmo:sub2lti:badInput',...
        ['I and J must be scalar integers or be\n' ...
        'arrays with the same number of integers!']);
end

% get first dimension if more than 1
nrows=nrows(1);

% subscripts should be in lower triangle
if(any(j>=i))
    error('seizmo:sub2lti:badInput','Indices outside lower triangle!');
end

% get lower triangle indices
k=cumsum(1:nrows-1);
lti=i(:)-k(j).'+nrows*(j(:)-1);

end
