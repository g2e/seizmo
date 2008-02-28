function [data]=chkhdr(data)
%CHKHDR    Check and update header field/values of seislab records
%
%    Description: Currently just updates header fields 'depmin', 'depmax',
%     'depmen', 'npts', and 'e' to reflect data records.  Field 'b' is
%     also updated for uneven records.
%
%    See also:  ch, gh, glgc, fixdelta

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% get header info
leven=glgc(data,'leven');
t=strcmp(leven,'true');
f=strcmp(leven,'false');
if(~all(t | f))
    error('sieslab:chkhdr:levenBad',...
        'logical field leven needs to be set'); 
end

% number of records
nrecs=length(data);

% preallocate
len=zeros(nrecs,1); depmax=zeros(nrecs,1);
depmin=zeros(nrecs,1); depmen=zeros(nrecs,1);

% loop over all records
for i=1:length(data)
    len(i)=size(data(i).x,1);
    depmax(i)=norm(max(data(i).x));
    depmin(i)=-norm(min(data(i).x));
    depmen(i)=norm(mean(data(i).x));
end

% work on uneven records
for i=find(f).'
    tlen=size(data(i).t,1);
    if(tlen~=len(i))
        error('seislab:chkhdr:dataSizeInconsistent',...
            ['number of samples in timing vector inconsistent with'...
            'data size for record %d'],i);
    end
    data(i)=ch(data(i),'b',data(i).t(1),'e',data(i).t(end));
end 

% work on even records
[b,delta]=gh(data(t),'b','delta'); 
data(t)=ch(data(t),'e',b+(len(t)-1).*delta);

% update some header fields
data=ch(data,'depmax',depmax,'npts',len,'depmin',depmin,'depmen',depmen);

end

