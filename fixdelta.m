function [data]=fixdelta(data,tol)
%FIXDELTA    Fix sample spacing for SEIZMO records
%
%    Usage:    data=fixdelta(data)
%
%    Description: DATA=FIXDELTA(DATA) fixes the sample spacing (DELTA
%     header field) of records in DATA to be the decimal equivalent of a
%     fraction of 2 small integers.  This is particularly useful for
%     upgrading the sample spacing from single to double precision as it
%     extends the precision.
%
%     DATA=FIXDELTA(DATA,TOL) allows specifying the maximum tolerance TOL
%     that the fraction of 2 small integers must match DELTA within.  So
%     a tolerance of 1e-4 requires the new sample rate to be within 1e-4 of
%     the old sample rate.
%
%    Notes:
%
%    Header changes: DELTA
%
%    Examples:
%     Force the double precision and update the delta field:
%      data=fixdelta(changeclass(data,'double'));
%
%    See also: rat, changeclass

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - minor doc update
%        Mar.  4, 2008 - minor doc update
%        Oct.  8, 2008 - doc update, add history
%        Nov. 16, 2008 - update for new name schema, doc and history update
%        Nov. 22, 2008 - doc update
%        Dec.  5, 2008 - tolerance option added
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% default tolerance
if(nargin==1 || ~isscalar(tol) || ~isnumeric(tol))
    [n,d]=rat(getheader(data,'delta'));
else
    [n,d]=rat(getheader(data,'delta'),tol);
end

% fix delta
data=changeheader(data,'delta',n./d);

end
