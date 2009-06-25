function [data]=keepph(data)
%KEEPPH    Returns phase component of spectral records
%
%    Usage:    data=keepph(data)
%
%    Description: KEEPPH(DATA) extracts the phase component of spectral
%     records in DATA.  Filetype is changed to General X vs Y.  This is
%     useful for operations that require only the phase component.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     Using keepph allows plotting the unwrapped phase using PLOT1:
%      plot1(seizmofun(keepph(data),@unwrap))
%
%    See also: keepam, keeprl, keepim, getspectralcmp, splitrecords, dft,
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
%     Last Updated June 25, 2009 at 16:55 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

data=getspectralcmp(data,'ph');

end
