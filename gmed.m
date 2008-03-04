function [med]=gmed(data)
%GMED    Returns median of each SAClab data record
%
%    Description: Returns the median of each SAClab data record.  Multi-
%     component records return a median for each component.  Output is
%     a cell array with one cell per record.
%
%    Usage: [medians]=gmed(data)
%
%    Examples:
%
%    See also: gnrm, rmean, subtr

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% find medians
med=cell(nrecs,1);
for i=1:nrecs
    med{i}=median(double(data(i).x));
end

end
