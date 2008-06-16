function [scale]=gnrm(data)
%GNRM    Returns normalizers for SAClab data records
%
%    Description: GNRM(DATA) returns the maximum amplitude of each SAClab 
%     data record in DATA assuming that for multi-component records the 
%     components are orthogonal to one another.  This is useful for 
%     normalizing multi-component data.  The maximum amplitude is just:
%                        _________________
%                       /
%       amp_max=max \  / cmp1^2+cmp2^2+...
%                    \/
%
%     Use header fields DEPMIN and DEPMAX to get the single largest value
%     in a record rather than a combined value.
%
%    Usage: scale=gnrm(data)
%
%    Examples:
%     Find the maximum amplitude of ground motion for a 3 component record:
%       amp_max=gnrm(data)
%
%    See also: mul, nrm, gmedian

%     Version History:
%        ????????????? - Initial Version
%        June 16, 2008 - Updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2008 at 03:25 GMT

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
