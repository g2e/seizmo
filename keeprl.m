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
%    See also: keepam, keepph, keepim, getspectralcmp, splitrecords, dft,
%              idft

%     Version History:
%        June 25, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 25, 2009 at 17:00 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

data=getspectralcmp(data,'rl');

end
