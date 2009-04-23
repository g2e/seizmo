function [data]=reverse(data)
%REVERSE    Reverse SEIZMO records
%
%    Usage:    data=reverse(data)
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
%    Examples:
%     See what this operation does:
%      plot1([data(1) reverse(data(1))])
%
%    See also:

%     Version History:
%        Apr.  9, 2008 - initial version
%        Nov. 22, 2008 - doc update, history fix, .dep
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 20:55 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% reverse records
for i=1:numel(data)
    data(i).dep=data(i).dep(end:-1:1,:);
end

end
