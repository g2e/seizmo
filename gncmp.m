function [ncmp]=gncmp(data)
%GNCMP    Returns the number of dependent components for each SAClab data record
%
%    Description: GNCMP(DATA) returns the number of dependent components
%     for each record in DATA.  This is basically a workaround for handling
%     SAClab versions that do not support multiple components (and thus not
%     having the NCMP field) at the same time as versions that do.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    ncmp=gncmp(data)
%
%    Examples:
%     Compare the number of components for a time series and spectral file:
%      gncmp([data(1) dft(data(1))])
%
%    See also: gh

%     Version History:
%        Oct.  7, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 03:35 GMT

% todo:

% input check
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'name','endian'))

% get ncmp via GH, avoiding warnings
warning('off','SAClab:gh:fieldInvalid')
ncmp=gh(data,'ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:gncmp:badNumCmp',...
        'Field NCMP must be a positive integer!')
end

end
