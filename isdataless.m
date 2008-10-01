function [lgc]=isdataless(data)
%ISDATALESS    True for dataless records in a SAClab structure
%
%    Description: ISDATALESS(DATA) returns a logical array equal in size to
%     DATA with elements set to true for the corresponding records in DATA
%     that are dataless and false otherwise.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: N/A
%
%    Usage:    lgc=isdataless(data)
%
%    Examples:
%     Skip dataless records in a for loop:
%      for i=find(~isdataless(data)); ...some content...; end
%
%    See also: hasmulpts, isseis

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

% check if dataless
lgc=false(size(data));
for i=1:numel(data)
    if(isempty(data(i).dep))
        lgc(i)=true;
    end
end

end