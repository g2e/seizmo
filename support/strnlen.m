function [A]=strnlen(A,n)
%STRNLEN    Pad/truncate char/cellstr array to n character columns
%
%    Usage: A=strnlen(A,len)
%
%    Description: STRNLEN(A,N) pads or truncates char/cellstr array A to
%     have N columns.  Padding is done with spaces.  Works recursively, so
%     will accept nested cellstr arrays.  Does not modify array type.
%
%    Notes:
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
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 28, 2009 - added testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 28, 2009 at 17:50 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

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
