function [data]=hilbrt(data)
%HILBRT    Return Hilbert transform of SEIZMO records
%
%    Description: HILBRT(DATA) calculates and returns the Hilbert transform
%     of SEIZMO records.  The Hilbert tranform is simply a -90 degree phase
%     shift.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=hilbrt(data)
%
%    Examples:
%     To do a positive 90 degree phase shift:
%      data=hilbrt(multiply(data,-1))
%
%    See also: envelope

%     Version History:
%        Jan. 30, 2008 - initial version
%        Feb. 28, 2008 - seischk support and class support
%        Mar.  4, 2008 - doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - minor doc update
%        July 17, 2008 - history update, doc update, now uses
%                        SEIZMO functions rather than Matlab's hilbert
%        Nov. 22, 2008 - doc update, update for new name schema
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 22:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin))

% transform, add phase shift, inverse transform
data=idft(subtract(dft(data),pi/2,2));

end
