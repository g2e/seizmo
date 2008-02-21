function [data]=distro(data,recmatrix,indices,store,npts)
%DISTRO    Distributes a record matrix back to SAClab structure
%
%    Description: Imports data records that are in a recmatrix into a
%     SAClab data structure.  Number of columns (records) in the record
%     matrix must match up with the number of records in the structure.
%     This is only for evenly spaced timeseries/xy files.  Column 1 goes
%     to data(1).x and so on.
%
%    Usage: [data]=distro(data,recmatrix,indices,store,npts)
%
%    See also: combo

% check input
error(nargchk(5,5,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% loop through records
for i=1:length(data)
    % retrieve record from matrix
    oclass=str2func(store{i});
    data(i).x=oclass(recmatrix(1:npts(i),i==indices));
end

% updata header
data=checkup(data);

end
