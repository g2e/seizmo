function [data]=doubleit(data)
%DOUBLEIT    Change in memory SAClab data storage to double precision
%
%    Description: Changes the data in memory to double precision.  This
%     does not affect the storage type or version for files written from 
%     the data (that requires changing the header version).
%
%    Usage: data=doubleit(data)
%
%    Examples:
%     Double the precision of records and fix the delta intervals
%      data=fixdelta(doubleit(data))
%
%    See also: combo, distro, fixdelta

% check nargin
error(nargchk(1,1,nargin))

% check data structure
if(~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% number of records
nrecs=length(data);

% retreive header info
[leven]=gh(data,'leven');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% add zeros and update header
for i=1:nrecs
    % logical index of header
    v=data(i).version==vers;
    
    % double header and data
    data(i).head=double(data(i).head);
    data(i).x=double(data(i).x);
    
    % double time
    if(leven(i)==h(v).false)
        data(i).t=double(data(i).t);
    end
end

end