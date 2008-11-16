function [median]=getmedian(data)
%GETMEDIAN    Returns median of each SEIZMO data record
%
%    Description: GETMEDIAN(DATA) returns the median of each SEIZMO data 
%     record.  Output is a cell array with one cell per record and multi-
%     component records return a vector of medians, one for each component.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    medians=getmedian(data)
%
%    Examples:
%     Remove median from records (cell2mat converts the cell array output
%     to a numeric matrix as long as all the records have the same number 
%     of components):
%      data=subtract(data,cell2mat(getmedian(data)));
%
%    See also: getnorm, removemean, subtract

%     Version History:
%        Feb. 12, 2008 - initial version
%        Feb. 28, 2008 - uses SEISCHK, doc update
%        Mar.  4, 2008 - minor doc update
%        June 15, 2008 - renamed from GMED to GMEDIAN, doc update
%        Nov. 16, 2008 - updated for new name schema (now GETMEDIAN),
%                        history fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 16, 2008 at 02:15 GMT

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% number of records
nrecs=numel(data);

% find medians
median=cell(nrecs,1);
for i=1:nrecs
    median{i}=median(double(data(i).dep));
end

end
