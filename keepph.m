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
%    See also: KEEPAM, KEEPRL, KEEPIM, GETSPECTRALCMP, SPLITRECORDS, DFT,
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

data=getspectralcmp(data,'ph');

end
