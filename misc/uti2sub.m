function [i,j]=uti2sub(nrows,uti)
%UTI2SUB    Square matrix upper triangle linear indices to subscripts
%
%    Usage:    [i,j]=uti2sub(nrows,idx)
%
%    Description:
%     [I,J]=UTI2SUB(NROWS,IDX) converts upper triangle indices IDX to
%     row/column subscripts I & J.  The number of rows in the matrix needs
%     to be given as NROWS.  Upper triangle indices increase similar to
%     linear indices, but the lower triangle and the diagonal are ignored
%     so that for a 5x5 matrix the indices proceed as follows:
%
%       \j  1  2  3  4  5
%      i \
%      1    -  1  2  4  7
%      2    -  -  3  5  8
%      3    -  -  -  6  9
%      4    -  -  -  - 10
%      5    -  -  -  -  -
%
%    Notes:
%
%    Examples:
%     % So you have a dissimilarity vector (in this case it corresponds
%     % to the upper triangle of a 400x400 dissimilarity matrix) and you
%     % wanted to know the indices of the two things being compared at
%     % index 8544 in the vector:
%     [thing1,thing2]=uti2sub(400,8544)
%
%    See also: LTI2SUB, SUB2UTI, SUB2LTI

%     Version History:
%        Sep.  8, 2009 - added documentation
%        Oct. 13, 2009 - added checks, updated documentation
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 22:35 GMT

% todo:

% checks
if(isempty(nrows) || ~isnumeric(nrows) ...
        || any(nrows~=fix(nrows)) || numel(nrows)>2)
    error('seizmo:uti2sub:badInput',...
        'NROWS must be a scalar number or 1x2 array of [NROWS NCOLS]!');
elseif(~isnumeric(uti) || any(uti~=fix(uti)))
    error('seizmo:uti2sub:badInput','IDX must be indices!');
end

% get first dimension if more than 1
nrows=nrows(1);

% get subscripts
uti=uti(:);
k=cumsum([0 1:nrows-1]);
[i,j]=min((uti(:,ones(1,nrows))>k(ones(numel(uti),1),:)).');
i=uti-k(j-1).';
j=j(:);

end
