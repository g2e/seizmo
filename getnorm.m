function [scale]=getnorm(data)
%GETNORM    Return normalizers for SEIZMO records
%
%    Usage:    scale=getnorm(data)
%
%    Description: GETNORM(DATA) returns the maximum amplitude of each
%     SEIZMO record in DATA assuming that for multi-component records
%     the components are orthogonal to one another.  This is useful for 
%     normalizing multi-component data.  The maximum amplitude is just:
%                        _________________
%                       /
%       amp_max=max \  / cmp1^2+cmp2^2+...
%                    \/
%
%     Use header fields DEPMIN and DEPMAX to get the single largest value
%     in a record rather than a combined value.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     Find the maximum amplitude of ground motion for a 3 component record:
%       amp_max=getnorm(data)
%
%    See also: multiply, normalize, getmedian

%     Version History:
%        Feb. 12, 2008 - initial version
%        Feb. 23, 2008 - bug fix
%        Feb. 28, 2008 - real bug fix
%        Mar.  4, 2008 - minor doc update
%        Apr. 17, 2008 - minor doc update
%        June 16, 2008 - doc update
%        Nov. 16, 2008 - renamed from GNRM to GETNORM
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 20:20 GMT

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% find normalizers
scale=nan(nrecs,1);
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    scale(i)=max(sqrt(sum(double(data(i).dep).^2,2)));
end

end
