function [data]=chkhdr(data)
%CHKHDR    Check header field/values of SAClab records
%
%    Description: Currently just updates header fields 'depmin', 'depmax',
%     'depmen', 'npts', and 'e'.
%
%    See also: fixdelta

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
len=zeros(nrecs,1); depmax=zeros(nrecs,1);
depmin=zeros(nrecs,1); depmen=zeros(nrecs,1);
for i=1:length(data)
    len(i)=size(data(i).x,1);
    depmax(i)=norm(max(data(i).x));
    depmin(i)=-norm(min(data(i).x));
    depmen(i)=norm(mean(data(i).x));
end

% get header and update
[b,delta]=gh(data,'b','delta'); e=b+(len-1).*delta;
data=ch(data,'depmax',depmax,'npts',len,'e',e,'depmin',depmin,'depmen',depmen);

end

