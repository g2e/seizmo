function [ncmp]=getncmp(data)
%GETNCMP    Return the number of dependent components for SEIZMO records
%
%    Description: GETNCMP(DATA) returns the number of dependent components
%     for each record in DATA.  This is basically a workaround for handling
%     filetypes that do not support multiple components (and thus not
%     having the NCMP field) at the same time as ones that do.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    ncmp=getncmp(data)
%
%    Examples:
%     Compare the number of components for a time series and spectral file:
%      getncmp([data(1) dft(data(1))])
%
%    See also: getheader, getenumid, getenumdesc, getlgc

%     Version History:
%        Oct.  7, 2008 - initial version
%        Nov. 16, 2008 - update for new name schema (now GETNCMP)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 03:35 GMT

% todo:

% input check
error(nargchk(1,1,nargin))

% get ncmp via GH, avoiding warnings
warning('off','seizmo:getheader:fieldInvalid')
ncmp=getheader(data,'ncmp');
warning('on','seizmo:getheader:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('seizmo:getncmp:badNumCmp',...
        'Field NCMP must be a positive integer!')
end

end
