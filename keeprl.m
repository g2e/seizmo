function [data]=keeprl(data)
%KEEPRL    Returns the real component of spectral records
%
%    Usage:    data=keeprl(data)
%
%    Description: KEEPRL(DATA) extracts the real component of spectral
%     records in DATA.  Filetype is changed to General X vs Y.  This is
%     useful for operations that require only the real component of the
%     spectra.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     Real vs imaginary:
%      p2([keeprl(data(1)) keepim(data(1))])
%
%    See also: KEEPAM, KEEPPH, KEEPIM, GETSPECTRALCMP, SPLITRECORDS, DFT,
%              IDFT

%     Version History:
%        June 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:25 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

data=getspectralcmp(data,'rl');

end
