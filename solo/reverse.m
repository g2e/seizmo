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
%    Header changes: NONE
%
%    Examples:
%     See what this operation does:
%      plot1([data(1) reverse(data(1))])
%
%    See also: MIRRORFLIP

%     Version History:
%        Apr.  9, 2008 - initial version
%        Nov. 22, 2008 - doc update, history fix, .dep
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Sep. 11, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% reverse records
for i=1:numel(data)
    data(i).dep=data(i).dep(end:-1:1,:);
end

end
