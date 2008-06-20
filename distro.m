function [data]=distro(data,recmatrix,indices,store,npts)
%DISTRO    Distributes a record matrix back to SAClab structure
%
%    Description: DISTRO(DATA,RECORDS,INDICES,STORE,NPTS) imports data in 
%     RECORDS (each component has its own column, INDICES gives the record
%     that each component belong to, STORE is the class the data should be
%     stored as in DATA, and NPTS gives the length of each record).  This
%     function is typically used in conjunction with COMBO to access
%     functions not in SAClab.  For creating a new SAClab data structure
%     see the command BSEIS.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: NONE
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, E, B
%
%    Usage: [data]=distro(data,recmatrix,indices,store,npts)
%
%    Examples:
%     Five sample running average:
%      [recmatrix,indices,store,npts]=combo(data);
%      recmatrix=filter(ones(5,1)/5,1,recmatrix);
%      data=distro(data,recmatrix,indices,store,npts);
%
%    See also: combo, bseis

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 03:30 GMT

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
