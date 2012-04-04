function [scale]=getnorm(data)
%GETNORM    Return normalizers for SEIZMO records
%
%    Usage:    scale=getnorm(data)
%
%    Description:
%     SCALE=GETNORM(DATA) returns the maximum amplitude of each SEIZMO
%     record in DATA assuming that for multi-component records the
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
%    Notes:
%
%    Examples:
%     % Find the maximum amplitude for a 3 component record:
%     amp_max=getnorm(data)
%
%    See also: MULTIPLY, NORMALIZE, GETVALUEFUN

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
%        Mar. 26, 2010 - doc update
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 15:05 GMT

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% find normalizers
scale=nan(nrecs,1);
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    scale(i)=max(sqrt(sum(double(data(i).dep).^2,2)));
end

end
