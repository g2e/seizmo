function [data]=reverse(data)
%REVERSE    Reverse SEIZMO records
%
%    Description: REVERSE(DATA) reverses the dependent components in
%     SEIZMO records so that the start and end of each is switched.
%
%    Notes:
%     - useful for flipping negative & positive frequencies
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    data=reverse(data)
%
%    Examples:
%     See what this operation does:
%      plot1([data(1) reverse(data(1))])
%
%    See also:

%     Version History:
%        Apr.  9, 2008 - initial version
%        Nov. 22, 2008 - doc update, history fix, .dep
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 09:20 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% reverse records
for i=1:numel(data)
    data(i).dep=data(i).dep(end:-1:1,:);
end

end
