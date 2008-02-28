function [data]=distro(data,recmatrix,indices,store,npts)
%DISTRO    Distributes a record matrix back to seislab structure
%
%    Description: Imports data records that are in a record matrix back
%     into a seislab data structure.
%
%    Usage: [data]=distro(data,recmatrix,indices,store,npts)
%
%    See also: combo

% check input
error(nargchk(5,5,nargin))

% check data structure
error(seischk(data))

% loop through records
for i=1:length(data)
    % retrieve record from matrix
    oclass=str2func(store{i});
    data(i).x=oclass(recmatrix(1:npts(i),i==indices));
end

% updata header
data=chkhdr(data);

end

