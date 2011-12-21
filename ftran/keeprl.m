function [data]=keeprl(data)
%KEEPRL    Returns the real component of spectral records
%
%    Usage:    data=keeprl(data)
%
%    Description:
%     KEEPRL(DATA) extracts the real component of spectral records in DATA.
%     Filetype is changed to General X vs Y.  This is useful for operations
%     that require only the real component of the spectra.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     % Real vs imaginary:
%     plot2([keeprl(data(1)) keepim(data(1))])
%
%    See also: KEEPAM, KEEPPH, KEEPIM, KEEPPW, GETSPECTRALCMP,
%              SPLITRECORDS, DFT, IDFT

%     Version History:
%        June 25, 2009 - initial version
%        Aug. 15, 2010 - nargchk fix
%        Dec. 21, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 21, 2011 at 20:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

data=getspectralcmp(data,'rl');

end
