function [median]=gmedian(data)
%GMEDIAN    Returns median of each SAClab data record
%
%    Description: GMEDIAN(DATA) returns the median of each SAClab data 
%     record.  Output is a cell array with one cell per record and multi-
%     component records return a vector of medians, one for each component.
%
%    Usage: [medians]=gmedian(data)
%
%    Examples:
%     Remove median from records (cell2mat converts the cell array output
%     to a numeric matrix as long as all the records have the same number 
%     of components):
%      data=sub(data,cell2mat(gmedian(data)));
%
%    See also: gnrm, rmean, sub

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Renamed and updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 04:05 GMT

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% find medians
median=cell(nrecs,1);
for i=1:nrecs
    median{i}=median(double(data(i).x));
end

end
