function [scale]=gnrm(data)
%GNRM    Returns normalizers for SAClab data records
%
%    Description: Returns the maximum amplitude of SAClab data records.
%
%    Usage: [scale]=gnrm(data)
%
%    See also: mul, nrm, gmed

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
    scale(i)=max(sqrt(sum((double(data(i).x).^2).',2)));
end

end

