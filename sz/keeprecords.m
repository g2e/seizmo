function [data]=keeprecords(data,idx)
%KEEPRECORDS    Keeps indicated records in SEIZMO data structure
%
%    Usage:    data=keeprecords(data,idx)
%
%    Description:
%     KEEPRECORDS(DATA,IDX) returns the SEIZMO structure with only the
%     records indicated by IDX.  IDX must be an array of linear indices or
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
%    See also: SELECTRECORDS, DELETERECORDS

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

data=data(idx);

end
