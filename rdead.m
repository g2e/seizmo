function [data]=rdead(data)
%RDEAD    Removes unvarying SAClab data records
%
%    Description: Removes dead records that can cause problems in analysis
%     and are not worth keeping.  Uses depmin/depmax so that records can be
%     eliminated before reading in the data.
%
%    Usage:  [data]=rdead(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: gh

% check input
error(nargchk(1,1,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% get depmax, depmin
[depmax,depmin]=gh(data,'depmax','depmin');

% remove dead records
data(depmax-depmin==0)=[];

end