function [data]=keepph(data)
%KEEPPH    Returns phase component of spectral records
%
%    Usage:    data=keepph(data)
%
%    Description:
%     KEEPPH(DATA) extracts the phase component of spectral records in
%     DATA.  Filetype is changed to General X vs Y.  This is useful for
%     operations that require only the phase component.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     % Using keepph allows plotting the unwrapped phase using PLOT1:
%     plot1(solofun(keepph(data),@unwrap))
%
%    See also: KEEPAM, KEEPRL, KEEPIM, KEEPPW, GETSPECTRALCMP,
%              SPLITRECORDS, DFT, IDFT

%     Version History:
%        June 25, 2009 - initial version
%        Aug. 15, 2010 - nargchk fix
%        Jan.  6, 2011 - fix name breakage in example
%        Dec. 21, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 21, 2011 at 20:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

data=getspectralcmp(data,'ph');

end
