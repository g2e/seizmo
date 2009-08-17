function [data]=envelope(data)
%ENVELOPE    Return envelopes of SEIZMO records
%
%    Usage:    data=envelope(data)
%
%    Description: ENVELOPE(DATA) returns the envelope (the complex
%     magnitude of a record's analytic signal) of the SEIZMO records.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Plot the envelopes against the data:
%      recordsection([data; envelope(data)])
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
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% sqrt(H(x)^2+x^2)
data=seizmofun(...
    addrecords(seizmofun(hilbrt(data),@(x)x.^2),...
    seizmofun(data,@(x)x.^2)),@(x)x.^0.5);

end
