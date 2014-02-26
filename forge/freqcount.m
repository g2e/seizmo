function [y,n,i]=freqcount(x)
%FREQCOUNT    Returns counts for values in X
%
%    Usage:    [y,n,i]=freqcount(x)
%
%    Description:
%     [Y,N,I]=FREQCOUNT(X) counts the number of times and the indices for
%     each unique value in X.  Y contains the unique values of X, N is the
%     number of occurances for each unique value and I is a cell array
%     containing the indices of the elements in X with those unique values.
%     Thus Y, N & I are all column vectors of equal size and the
%     corresponding elements give the details for a particular unique value
%     in X.
%
%    Notes:
%     - NaNs are handled as equal and are placed at the end of the outputs
%       if there are any.
%     - Also works with char & cellstr arrays.  Logical arrays are
%       supported but there is no point to using FREQCOUNT for that type.
%
%    Examples:
%     % Demonstrate on a Hankel matrix:
%     x=hankel(1:10);
%     [y,n,i]=freqcount(x)
%     x(i{3}) % directly index the 2s of x
%
%     % Character counts for this function:
%     [y,n,i]=freqcount(readtxt(which('freqcount')));
%     [n,j]=sort(n);
%     tmp=[num2cell(n) cellstr(y(j))]';
%     sprintf('%d %s\n',tmp{:})
%
%     % Word frequency for this function:
%     [y,n,i]=freqcount(getwords(readtxt(which('freqcount'))));
%     [n,j]=sort(n);
%     tmp=[num2cell(n) y(j)]';
%     sprintf('%d %s\n',tmp{:})
%
%    See also: MODE, HISTC

%     Version History:
%        Feb. 20, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 11:15 GMT

% todo:

% require an input
error(nargchk(1,1,nargin));

% convert to column vector
x=x(:);

% numeric inputs avoid unique for speed
if(isnumeric(x))
    % handle nans beforehand
    nans=isnan(x);
    nn=find(~nans);
    
    % avoid unique b/c it is slow
    [xs,i]=sort(x(~isnan(x)));
    y=xs([find(diff(xs));end]);
    n=histc(xs,y);
    i=mat2cell(nn(i),n);
    
    % account for nans
    if(any(nans))
        y(end+1)=nan;
        i{end+1}=find(nans);
        n(end+1)=numel(i{end});
    end
else % cell string
    [y,i,j]=unique(x);
    n=histc(j,1:numel(y));
    lind=1:numel(x);
    [i,i]=sort(j);
    i=mat2cell(lind(i)',n);
end

end
