function [data]=keepam(data)
%KEEPAM    Returns amplitudes component of spectral records
%
%    Usage:    data=keepam(data)
%
%    Description: KEEPAM(DATA) extracts the amplitude component of spectral
%     records in DATA.  Filetype is changed to General X vs Y.  This is
%     useful for operations that require only spectral amplitudes.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     Using keepam allows plotting amplitudes using PLOT1:
%      plot1(keepam(data))
%
%    See also: keepph, keeprl, keepim, getspectralcmp, splitrecords, dft,
%              idft

%     Version History:
%        June 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:25 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

data=getspectralcmp(data,'am');

end
