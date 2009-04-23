function [median]=getmedian(data)
%GETMEDIAN    Returns median of each SEIZMO record
%
%    Usage:    medians=getmedian(data)
%
%    Description: GETMEDIAN(DATA) returns the median value in the dependent
%     component(s) for each SEIZMO record.  Multi-component records are
%     just combined and treated as a single one to find the median.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     Remove median from records:
%      data=subtract(data,getmedian(data));
%
%    See also: getnorm, removemean, subtract

%     Version History:
%        Feb. 12, 2008 - initial version
%        Feb. 28, 2008 - uses SEISCHK, doc update
%        Mar.  4, 2008 - minor doc update
%        June 15, 2008 - renamed from GMED to GMEDIAN, doc update
%        Nov. 16, 2008 - updated for new name schema (now GETMEDIAN),
%                        history fix
%        Nov. 23, 2008 - output changed to just one value per record not
%                        one value per component
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 20:15 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% find medians
median=nan(nrecs,1);
for i=1:nrecs
    median(i)=median(double(data(i).dep(:)));
end

end
