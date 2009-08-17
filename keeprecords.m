function [data]=keeprecords(data,idx)
%KEEPRECORDS    Keeps indicated records in SEIZMO data structure
%
%    Usage:    data=keeprecords(data,idx)
%
%    Description: KEEPRECORDS(DATA,IDX) returns the SEIZMO structure
%     with only the records indicated by IDX.  IDX must be an array of
%     linear indices or a logical index array.  This could be useful for
%     inline record deletion.
%
%    Notes:
%
%    Examples:
%     These do the exact same thing:
%      data=keeprecords(data,[1 3 6:end]);
%      data=data([1 3 6:end]);
%      data=deleterecords(data,[2 5]);
%      data([2 5])=[];
%
%    See also: selectrecords, deleterecords

%     Version History:
%        June 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:25 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

data=data(idx);

end
