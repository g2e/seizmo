function [i,j]=uti2sub(len,uti)
%UTI2SUB    Square matrix upper triangle linear indices to subscripts
%
%    Usage:    [i,j]=uti2sub(nrows,idx)
%
%    Description: [I,J]=uti2sub(NROWS,IDX) converts upper triangle indices
%     IDX to row/column subscripts I & J.  The number of rows in the matrix
%     needs to be given as NROWS.  Upper triangle indices increase similar
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
%     - No input checks are done!
%
%    Examples:
%     Say you have a dissimilarity vector (in this case, say it corresponds
%     to the upper triangle of a 400x400 dissimilarity matrix) and you
%     wanted to know the two things that are being compared at value 8544
%     in the vector:
%      [thing1,thing2]=uti2sub(400,8544)
%
%    See also: LTI2SUB, SUB2UTI, SUB2LTI

%     Version History:
%        Sep.  8, 2009 - added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 04:50 GMT

% todo:

len=len(1);
uti=uti(:);
k=cumsum([0 1:len-1]);
[i,j]=min((uti(:,ones(1,len))>k(ones(length(uti),1),:)).');
i=uti-k(j-1).'; j=j(:);

end
