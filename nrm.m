function [data,scale]=nrm(data)
%NRM    Normalizes SAClab data records
%
%    Description: Normalizes the amplitudes of SAClab data records to the
%     range of -1 to 1.  Second output contains the normalization factors.
%
%    Usage: [data,scale]=nrm(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: mul

% check nargin
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

% normalize data
scale=ones(nrecs,1);
for i=1:nrecs
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    scale(i)=max(sqrt(sum(((data(i).x).^2).')));
    data(i).x=oclass(data(i).x/scale(i));
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
