function [data]=keepim(data)
%KEEPIM    Returns the imaginary component of spectral records
%
%    Usage:    data=keepim(data)
%
%    Description:
%     KEEPIM(DATA) extracts the imaginary component of spectral records in
%     DATA.  Filetype is changed to General X vs Y.  This is useful for
%     operations that require only the imaginary component of the spectra.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     % Real vs imaginary:
%     plot2([keeprl(data(1)) keepim(data(1))])
%
%    See also: KEEPAM, KEEPPH, KEEPRL, KEEPPW, GETSPECTRALCMP,
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

data=getspectralcmp(data,'im');

end
