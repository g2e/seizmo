function [data]=rdead(data)
%RDEAD    Removes unvarying seislab data records
%
%    Description: Removes dead records that can cause problems in analysis
%     and are not worth keeping.  Uses depmin/depmax so that records can be
%     eliminated before reading in the data.
%
%    Usage:  [data]=rdead(data)
%
%    See also: gh

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% get depmax, depmin
[depmax,depmin]=gh(data,'depmax','depmin');

% remove dead records
data(depmax-depmin==0)=[];

end
