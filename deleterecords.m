function [data]=deleterecords(data,idx)
%DELETERECORDS    Deletes indicated records from SEIZMO data structure
%
%    Usage:    data=deleterecords(data,idx)
%
%    Description: DELETERECORDS(DATA,IDX) returns the SEIZMO structure
%     with the records indicated by IDX removed.  IDX must be an array of
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
%    See also: selectrecords, keeprecords

%     Version History:
%        June 25, 2009 - initial version
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
%     Last Updated June 25, 2009 at 09:15 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

data(idx)=[];

end
