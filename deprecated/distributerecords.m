function [data]=distro(data,recmatrix,indices,store,npts)
%DISTRO    Distributes a record matrix back to SAClab structure
%
%    Description: DISTRO(DATA,RECORDS,INDICES,STORE,NPTS) imports the data 
%     in RECORDS (a numeric array with each column as a separate component
%     or record) into SAClab structure DATA and outputs the updated 
%     structure.  INDICES gives the record index that each column 
%     corresponds to (if multiple columns belong to a record they will be 
%     ordered the same as they are in RECORDS).  STORE gives the data 
%     storage class for each record (not for each column/component).  NPTS 
%     gives the number of points in each record (not for each 
%     column/component).  This function is typically used in conjunction 
%     with COMBO to access functions not in SAClab.  For creating a new 
%     SAClab data structure see the command BSEIS.
%
%    Notes:
%     - INDICES should have as many elements as RECORDS has columns
%     - STORE should be a cell string array with as many elements as DATA
%     - NPTS should have as many elements as DATA
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
%        Feb. 16, 2008 - initial version
%        Feb. 21, 2008 - minor documentation update
%        Feb. 23, 2008 - minor documentation update
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor documentation update
%        June 15, 2008 - documentation update
%        June 29, 2008 - documentation update, .dep rather than .x
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 29, 2008 at 06:15 GMT

% todo:
% - option to import independent component too

% check input
error(nargchk(5,5,nargin))

% check data structure
error(seischk(data))

% loop through records
for i=1:numel(data)
    % retrieve record from matrix
    oclass=str2func(store{i});
    data(i).dep=oclass(recmatrix(1:npts(i),i==indices));
end

% updata header
data=chkhdr(data);

end
