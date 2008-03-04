function [data]=doubleit(data)
%DOUBLEIT    Change in memory SAClab data storage to double precision
%
%    Description: Changes the SAClab data in memory to double precision.
%     This does not affect the storage type or version for files written
%     from the data (that requires changing the header version).
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
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% retreive header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))

% add zeros and update header
for i=1:nrecs
    % double header and data
    data(i).head=double(data(i).head);
    data(i).x=double(data(i).x);
    
    % double time
    if(strcmp(leven(i),'false'))
        data(i).t=double(data(i).t);
    end
end

end
