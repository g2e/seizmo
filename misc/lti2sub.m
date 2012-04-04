function [i,j]=lti2sub(nrows,lti)
%LTI2SUB    Square matrix lower triangle linear indices to subscripts
%
%    Usage:    [i,j]=lti2sub(nrows,idx)
%
%    Description:
%     [I,J]=LTI2SUB(NROWS,IDX) converts lower triangle indices IDX to
%     row/column subscripts I & J.  The number of rows in the matrix needs
%     to be given as NROWS.  Lower triangle indices increase similar to
%     linear indices, but the upper triangle and the diagonal are ignored
%     so that for a 5x5 matrix the indices proceed as follows:
%
%       \j  1  2  3  4  5
%      i \
%      1    -  -  -  -  -
%      2    1  -  -  -  -
%      3    2  5  -  -  -
%      4    3  6  8  -  -
%      5    4  7  9 10  -
%
%    Notes:
%
%    Examples:
%     % So you have a dissimilarity vector (in this case it corresponds
%     % to the lower triangle of a 400x400 dissimilarity matrix) and you
%     % wanted to know the indices of the two things being compared at
%     % index 8544 in the vector:
%     [thing1,thing2]=lti2sub(400,8544)
%
%    See also: UTI2SUB, SUB2UTI, SUB2LTI

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
    error('seizmo:lti2sub:badInput',...
        'NROWS must be a scalar number or 1x2 array of [NROWS NCOLS]!');
elseif(~isnumeric(lti) || any(lti~=fix(lti)))
    error('seizmo:lti2sub:badInput','IDX must be indices!');
end

% get first dimension if more than 1
nrows=nrows(1);

% subscripts should be in lower triangle
if(any(lti>sum(1:nrows-1)))
    error('seizmo:lti2sub:badInput','Indices outside lower triangle!');
end

% get subscripts
lti=lti(:).';
k=cumsum([0 nrows-1:-1:1]).';
[i,j]=min((lti(ones(1,nrows),:)>k(:,ones(numel(lti),1))));
j=j(:)-1;
k=cumsum(1:nrows-1).';
i=lti(:)+k(j)-nrows*(j-1);

end
