function [med]=gmed(data)
%GMED    Returns medians of SAClab data records
%
%    Description: Returns the median of each SAClab data record.  Multi-
%     component records return a median for each component.  Output is
%     a cell array with one cell per record.
%
%    Usage: [medians]=gmed(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: gnrm, rmean, add

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
med=cell(nrecs,1);
for i=1:nrecs
    med{i}=median(data(i).x);
end

end