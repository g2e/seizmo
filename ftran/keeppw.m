function [data]=keeppw(data)
%KEEPPW    Returns the power component of spectral records
%
%    Usage:    data=keeppw(data)
%
%    Description:
%     KEEPPW(DATA) extracts the power component of spectral records in
%     DATA.  Filetype is changed to General X vs Y.  This is useful for
%     operations that require the power component of the spectra.  Units
%     are decibels (dBs).
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     % Simple power plot:
%     plot2(keeppw(data),'xscale','log');
%
%    See also: KEEPAM, KEEPPH, KEEPIM, KEEPRL, GETSPECTRALCMP,
%              SPLITRECORDS, DFT, IDFT

%     Version History:
%        Dec. 21, 2011 - initial version
%        Feb.  6, 2012 - power spectra is now in dBs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  6, 2012 at 20:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

data=getspectralcmp(data,'pw');

end
