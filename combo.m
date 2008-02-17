function [recmatrix,indices,store,npts]=combo(data)
%COMBO    Combines data records into a single numeric matrix
%
%    Description: Combines data records into one matrix (records are in
%     separate columns).  Records are padded with zeros to have the same
%     number of points.  This is for vectorizing heavy computations in the
%     short term.  Look into using multi-component files for better 
%     performance on large sets of records.
%
%    See also: distro

% check input
error(nargchk(1,1,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

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

% preallocate record matrix
col2=cumsum(ncol); col=[1; col2+1];
recmatrix=zeros(max(npts),col2(end));
indices=zeros(1,col2(end));

% loop through records
for i=1:nrecs
    indices(col(i):col2(i))=i;
    recmatrix(1:npts(i),col(i):col2(i))=data(i).x;
end

end
