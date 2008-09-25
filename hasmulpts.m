function [lgc]=hasmulpts(data)
%hasmulpts    True for records in a SAClab structure with multiple points
%
%    Description: HASMULPTS(DATA) returns a logical array equal in size to
%     DATA with elements set to true for the corresponding records in DATA
%     that have 2 or more data points and false otherwise.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: DATA must have 'dep' field
%
%    Header changes: N/A
%
%    Usage:    lgc=hasmulpts(data)
%
%    Examples:
%     Skip dataless and 1point records in a for loop:
%      for i=find(hasmulpts(data)); ...some content...; end
%
%    See also: isdataless, isseis

%     Version History:
%        Sep. 25, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2008 at 04:15 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% check if 2+ points
lgc=false(size(data));
for i=1:numel(data)
    if(size(data(i).dep,1)>1)
        lgc(i)=true;
    end
end

end