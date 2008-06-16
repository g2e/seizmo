function [data,scale]=nrm(data)
%NRM    Normalizes SAClab data records
%
%    Description: Normalizes the amplitudes of SAClab data records to the
%     range of -1 to 1.  Second output contains the normalization factors.
%     Components are assumed to be orthogonal when finding the
%     normalization.  Use header fields DEPMIN and DEPMAX to get the
%     single largest value.
%
%    Usage: [data,scale]=nrm(data)
%
%    Examples:
%
%    See also: mul, gnrm

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% normalize data
scale=ones(nrecs,1);
for i=1:nrecs
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    scale(i)=max(sqrt(sum((data(i).x).^2,2)));
    data(i).x=oclass(data(i).x/scale(i));
    data(i)=ch(data(i),'depmax',max(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmen',mean(data(i).x(:)));
end

end
