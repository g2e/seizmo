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
%    Usage: [recmatrix,indices,store,npts]=combo(data)
%
%    Examples:
%     Get the interquartile range of records and assign to a header field:
%      recs=combo(data);
%      data=ch(data,'user3',iqr(recs));
%
%    See also: distro

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 02:50 GMT

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% loop through records
store=cell(nrecs,1);
npts=zeros(nrecs,1);
ncol=zeros(nrecs,1);
for i=1:nrecs
    store{i}=class(data(i).x);
    [npts(i),ncol(i)]=size(data(i).x);
end

% preallocate record matrix (doubles)
col2=cumsum(ncol); col=[1; col2+1];
recmatrix=zeros(max(npts),col2(end));
indices=zeros(1,col2(end));

% loop through records
for i=1:nrecs
    indices(col(i):col2(i))=i;
    recmatrix(1:npts(i),col(i):col2(i))=data(i).x;
end

end
