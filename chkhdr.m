function [data]=chkhdr(data,varargin)
%CHKHDR    Check and update header field/values of SAClab records
%
%    Description: Currently just updates header fields 'depmin', 'depmax',
%     'depmen', 'npts', and 'e' to reflect data records.  Field 'b' is
%     also updated for uneven records.  Additional header changes can be
%     done by adding them as additional arguments.
%
%    Usage: data=chkhdr(data)
%           data=chkhdr(data,'field1',values1,'field2',values2,...)
%
%    Examples:
%     Check some basic header fields and also modify the t0 field with 
%     some arrivals:
%      data=chkhdr(data,'t0',P_arrivaltimes)
%
%    See also:  ch, gh, glgc, fixdelta

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% get header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
tru=strcmp(leven,'true');
fals=strcmp(leven,'false');

% number of records
nrecs=length(data);

% preallocate
len=zeros(nrecs,1); 
depmax=len; depmin=len; depmen=len;

% loop over all records, get npts, dep(min,max,men)
for i=1:nrecs
    len(i)=size(data(i).x,1);
    depmax(i)=max(data(i).x(:));
    depmin(i)=min(data(i).x(:));
    depmen(i)=mean(data(i).x(:));
end

% work only on uneven records, updating b and e fields
for i=find(fals).'
    tlen=size(data(i).t,1);
    if(tlen~=len(i))
        error('SAClab:chkhdr:dataSizeInconsistent',...
            ['number of samples in timing vector inconsistent with'...
            'data size for record %d'],i);
    end
    data(i)=ch(data(i),'b',data(i).t(1),'e',data(i).t(end));
end 

% work only on even records, updating e field
[b,delta]=gh(data(tru),'b','delta'); 
data(tru)=ch(data(tru),'e',b+(len(tru)-1).*delta);

% update some header fields for all records
data=ch(data,'depmax',depmax,'npts',len,...
    'depmin',depmin,'depmen',depmen,varargin{:});

end
