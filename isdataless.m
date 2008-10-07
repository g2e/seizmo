function [lgc]=isdataless(data)
%ISDATALESS    True for records in a SAClab structure without data
%
%    Description: ISDATALESS(DATA) returns a logical array equal in size to
%     DATA with elements set to TRUE for the corresponding records in DATA
%     that are dataless and false otherwise.  This is intended for checking
%     for records that have 0 points or have no data read in yet.
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
%     To find 0pt records:
%      ~gh(data,'npts')
%
%    See also: isseis

%     Version History:
%        Sep. 25, 2008 - initial version
%        Oct.  3, 2008 - fix .dep-less case, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  3, 2008 at 06:10 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% catch the case where no records have data read in
lgc=true(size(data));
if(isfield(data,'dep'))
    % check for data
    for i=1:numel(data)
        if(~isempty(data(i).dep))
            lgc(i)=false;
        end
    end
end

end