function [recmatrix,indices,store,npts]=combo(data)
%COMBO    Combines SAClab data records into a single numeric record matrix
%
%    Description: [RECORDS,INDICES,STORE,NPTS]=COMBO(DATA) combines SAClab
%     data records stored in DATA into one matrix RECORDS (records are in
%     separate columns).  INDICES indicates which columns belong to which 
%     records in DATA (useful when dealing with records with multiple 
%     components).  STORE is a cell array of data class strings, one for 
%     each record.  NPTS gives the original npts in each record (records 
%     are padded with zeros when combined).  This function is useful for 
%     providing easy access to functions not in SAClab.  Use DISTRO to
%     redistribute the records back into DATA.
%    
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: NONE
%
%    Header changes: NONE
%
%    Usage: [recmatrix,indices,store,npts]=combo(data)
%
%    Examples:
%     Get the interquartile range of records and assign to a header field:
%      recs=combo(data);
%      data=ch(data,'user3',iqr(recs));
%
%    See also: distro

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 21, 2008 - minor documentation update
%        Feb. 23, 2008 - minor documentation update
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor documentation update
%        June 15, 2008 - documentation update
%        June 20, 2008 - minor documentaiton update
%        June 29, 2008 - documentation update, .dep rather than .x,
%                        dataless support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 29, 2008 at 05:50 GMT

% todo:
% - option to extract independent component too

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% number of records
nrecs=numel(data);

% loop through records
store=cell(nrecs,1);
npts=zeros(nrecs,1);
ncol=zeros(nrecs,1);
for i=1:nrecs
    store{i}=class(data(i).dep);
    [npts(i),ncol(i)]=size(data(i).dep);
    if(npts(i)==0); ncol(i)=0; end
end

% preallocate record matrix (doubles)
col2=cumsum(ncol); col=[1; col2+1];
recmatrix=zeros(max(npts),col2(end));
indices=zeros(1,col2(end));

% loop through records
for i=1:nrecs
    indices(col(i):col2(i))=i;
    recmatrix(1:npts(i),col(i):col2(i))=data(i).dep;
end

end
