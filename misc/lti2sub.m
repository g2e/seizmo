function [i,j]=lti2sub(len,lti)
%LTI2SUB    Square matrix lower triangle linear indices to subscripts
%
%    Usage:    [i,j]=lti2sub(nrows,idx)
%
%    Description: [I,J]=lti2sub(NROWS,IDX) converts lower triangle indices
%     IDX to row/column subscripts I & J.  The number of rows in the matrix
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
%     - No input checks are done!
%
%    Examples:
%     Say you have a dissimilarity vector (in this case, say it corresponds
%     to the lower triangle of a 400x400 dissimilarity matrix) and you
%     wanted to know the two things that are being compared at value 8544
%     in the vector:
%      [thing1,thing2]=lti2sub(400,8544)
%
%    See also: UTI2SUB, SUB2UTI, SUB2LTI

%     Version History:
%        Sep.  8, 2009 - added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 04:50 GMT

% todo:

len=len(1);
if(any(idx>sum(1:len-1)))
    error('seizmo:lti2sub:badInput','Indices out of range!');
end
lti=lti(:);
k=cumsum([0 len-1:-1:1]);
[i,j]=min((lti(:,ones(1,len))>k(ones(length(lti),1),:)).');
j=j(:)-1; k=cumsum(1:len-1); i=lti+k(j).'-len*(j-1);

end
