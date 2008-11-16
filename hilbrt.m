function [data]=hilbrt(data)
%HILBRT    Return Hilbert transform of SEIZMO data records
%
%    Description: Calculates and returns the Hilbert transform of SEIZMO
%     data records.  The Hilbert tranform is simply a -90 degree phase 
%     shift.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: Evenly spaced; Time Series or General X vs Y
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=hilbrt(data)
%
%    Examples:
%     To do a positive 90 degree phase shift:
%      data=hilbert(mul(data,-1))
%
%    See also: envelope, mul

%     Version History:
%        Jan. 30, 2008 - initial version
%        Feb. 28, 2008 - seischk support and class support
%        Mar.  4, 2008 - doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - minor doc update
%        July 17, 2008 - history update, doc update, now uses
%                        SEIZMO functions rather than Matlab's hilbert
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2008 at 17:10 GMT

% todo:
%

% check nargin
error(nargchk(1,1,nargin))

% transform, add phase shift, inverse transform
data=idft(sub(dft(data),pi/2,2));

end
