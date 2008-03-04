function [scale]=gnrm(data)
%GNRM    Returns normalizers for SAClab data records
%
%    Description: Returns the maximum amplitude of SAClab data records.
%
%    Usage: [scale]=gnrm(data)
%
%    Examples:
%
%    See also: mul, nrm, gmed

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% find normalizers
scale=ones(nrecs,1);
for i=1:nrecs
    scale(i)=max(sqrt(sum(double(data(i).x).^2,2)));
end

end
