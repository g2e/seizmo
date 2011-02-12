function [data]=deleterecords(data,idx)
%DELETERECORDS    Deletes indicated records from SEIZMO data structure
%
%    Usage:    data=deleterecords(data,idx)
%
%    Description:
%     DELETERECORDS(DATA,IDX) returns the SEIZMO structure with the records
%     indicated by IDX removed.  IDX must be an array of linear indices or
%     a logical index array.  This could be useful for inline record
%     deletion.
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     % These all do the exact same thing:
%     data=keeprecords(data,[1 3 6:numel(data)]);
%     data=data([1 3 6:end]);
%     data=deleterecords(data,[2 5]);
%     data([2 5])=[];
%
%    See also: SELECTRECORDS, KEEPRECORDS

%     Version History:
%        June 25, 2009 - initial version
%        Oct. 21, 2009 - updated example
%        Feb. 11, 2011 - mass nargchk fix, minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

data(idx)=[];

end
