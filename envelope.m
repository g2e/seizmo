function [data]=envelope(data)
%ENVELOPE    Return envelopes of SEIZMO data records
%
%    Description: Returns the envelopes of the SEIZMO data records (complex
%     magnitude of a record's analytic signal).
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: Evenly spaced; Time Series or General X vs Y
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=envelope(data)
%
%    Examples:
%     Plot the envelopes against the data:
%      recsec([data; envelope(data)])
%
%    See also: hilbrt

%     Version History:
%        Jan. 30, 2008 - initial version
%        Feb. 23, 2008 - seischk support and class support
%        Feb. 28, 2008 - support for changing fft zeropadding
%        Mar.  4, 2008 - doc update
%        May  12, 2998 - dep* fix
%        July 17, 2008 - history update, doc update, now uses
%                        SEIZMO functions rather than Matlab's hilbert
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2008 at 17:10 GMT

% todo:
%

% check nargin
error(nargchk(1,1,nargin))

% sqrt(H(x)^2+x^2)
data=seisfun(...
    addf(seisfun(hilbrt(data),@(x)x.^2),seisfun(data,@(x)x.^2)),...
    @(x)x.^0.5);

end
