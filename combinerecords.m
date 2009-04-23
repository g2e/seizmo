function [dep,idx1,ind,idx2,store,npts]=combinerecords(data)
%COMBINERECORDS    Combines SEIZMO data records into a single numeric record matrix
%
%    Usage:    [dep,idx1,ind,idx2,store,npts]=combo(data)
%
%    Description: [DEP,IDX1,IND,IDX2,STORE,NPTS]=COMBINERECORDS(DATA)
%     combines SEIZMO records stored in DATA into numeric arrays DEP & IND 
%     (components are in separate columns).  IDX1 & IDX2 indicate which
%     columns belong to which records in DATA.  STORE is a cell array of
%     data class strings, one for each record.  NPTS gives the original
%     npts in each record (records are padded with zeros when combined).
%     This function is useful for providing easy access to functions not
%     in the SEIZMO toolbox.  Use DISTRIBUTERECORDS to redistribute the
%     records back into DATA if necessary.
%    
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: see CHECKHEADER
%
%    Examples:
%     Get the interquartile range of records and assign to a header field:
%      dep=combinerecords(data);
%      data=changeheader(data,'user3',iqr(dep));
%
%    See also: distributerecords, bseizmo

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 21, 2008 - minor doc update
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        June 15, 2008 - doc update
%        June 20, 2008 - minor doc update
%        June 29, 2008 - doc update, .dep rather than .x,
%                        dataless support
%        Nov. 22, 2008 - update for new name schema (now COMBINERECORDS)
%                        now outputs independent component
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 20:00 GMT

% todo:

% check input
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% get header fields
[b,npts,delta]=getheader(data,'b','npts','delta');

% number of records
nrecs=numel(data);

% loop through records
store=cell(nrecs,1);
npts=zeros(nrecs,1);
ncol=zeros(nrecs,1);
for i=1:nrecs
    store{i}=class(data(i).dep);
    [npts(i),ncol(i)]=size(data(i).dep);
    if(npts(i)*ncol(i)==0); npts(i)=0; ncol(i)=0; end
end
leven=true(nrecs,1);
if(isfield(data,'ind'))
    for i=1:nrecs
        if(~isempty(data(i).ind))
            leven(i)=false;
        end
    end
end

% preallocate dep matrix (as double precision)
col2=cumsum(ncol); col=[1; col2+1];
dep=zeros(max(npts),col2(end));
idx1=zeros(1,col2(end));

% preallocate ind matrix
ind=zeros(max(npts),nrecs);
idx2=1:nrecs;

% loop through records
for i=1:nrecs
    idx1(col(i):col2(i))=i;
    dep(1:npts(i),col(i):col2(i))=data(i).dep;
    if(leven(i))
        ind(1:npts(i),i)=b(i)+(0:npts(i)-1)*delta(i);
    else
        ind(1:npts(i),i)=data(i).ind;
    end
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
