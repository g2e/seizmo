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
%    See also: KEEPPH, KEEPRL, KEEPIM, GETSPECTRALCMP, SPLITRECORDS, DFT,
%              IDFT

%     Version History:
%        June 25, 2009 - initial version
%        Aug. 15, 2010 - nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 15, 2010 at 20:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

data=getspectralcmp(data,'am');

end
