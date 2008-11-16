function [A]=strnlen(A,n)
%STRNLEN    Pad/truncate char/cellstr array to n character columns
%
%    Description: STRNLEN(A,N) pads or truncates char/cellstr array A to
%     have N columns.  Padding is done with spaces.  Works recursively, so
%     will accept nested cellstr arrays.  Does not modify array type.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage: A=strnlen(A,len)
%
%    Examples:
%     Require elements in your array to be strings of length 8:
%      A=strnlen(A,8);
%
%    See also: strtrim, deblank, strjust

%     Version History:
%        Mar.  7, 2008 - initial version
%        Apr. 18, 2008 - fixed recursion
%        Oct. 26, 2008 - code cleaning, comment and doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 26, 2008 at 14:45 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin))

% check n
if(~isnumeric(n) || ~isscalar(n) || fix(n)~=n)
    error('seizmo:strnlen:badN','N must be an integer!');
end

% find elements to trucate/pad
if(ischar(A))
    d=size(A);
    if(d(2)>=n)
        % truncate
        A=A(:,1:n);
    else
        % pad with spaces
        A=[A ones(d(1),n-d(2))*32];
    end
elseif(iscellstr(A) || iscell(A))
    % recurse for each cell
    for i=1:numel(A)
        A{i}=strnlen(A{i},n);
    end
else
    error('seizmo:strnlen:badA',...
        'A must have only char elements!');
end

end
