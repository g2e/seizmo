function [data]=slidingam(data,nsamples,varargin)
%SLIDINGAM    Returns sliding-window absolute-mean of SEIZMO records
%
%    Description: SLIDINGAM(DATA,N) applies a centered sliding-window 
%     absolute-mean of 2N+1 samples to the dependent component(s) of 
%     SEIZMO data records in DATA.  N can be a scalar (each record has the
%     same window size) or a vector (define each record's window size 
%     separately).  Sliding windows extending outside the record are
%     truncated (look at 'EDGE' option to change this).
%
%     SLIDINGAM(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the position
%     of the sliding window relative to the reference data point.  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1, while for TRAIL or LEAD the window size
%     is N.  Default position is CENTER.
%     
%     SLIDINGAM(...,'OFFSET',OFFSET) sets the offset of the sliding window
%     from to the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of 
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     OFFSET may be a vector of offsets specifying each record's offset.
%     
%     SLIDINGAM(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the edge
%     cases.  TRUNCATE eliminates points in the sliding-window that do not 
%     reference a datapoint (for instance if the window extends before or 
%     after the data, that portion of the window will be truncated).  PAD 
%     adds zeros to the data so that all the points of the sliding-window 
%     always reference some value.  Default setting is TRUNCATE.
%
%     SLIDINGAM(...,'DIM',N) specifies an alternative dimension to slide
%     across rather than the default 1 (slides down the component - 2 would
%     slide across the components).
%
%     SLIDINGAM(...,'CUSTOM',WINDOW) allows a custom sliding window
%     average.  This might be useful for a Gaussian average or similar.
%     WINDOW must be formatted as [index; weight] where index is relative
%     to the reference data point and weight does not include the averaging
%     divisor (this will be automatically computed).  An example WINDOW:
%        [ -3  -2  -1    0   1  2    3;
%         0.1   1  10  100  10  1  0.1]
%     gives a severe weighting to the center point (reference data point).
%
%    Notes:
%     - Centered windows are of length 2N+1, while the others are just N
%     - SLIDINGAM is faster than SLIDEFUN because it uses SLIDINGAVG
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=slidingam(data,n)
%              data=slidingam(...,'position','center'|'trail'|'lead')
%              data=slidingam(...,'offset',offset)
%              data=slidingam(...,'edge','truncate'|'pad')
%              data=slidingam(...,'dim',n)
%              data=slidingam(...,'custom',window)
%
%    Examples:
%      Compare an envelope and a 21-sample sliding-window absolute-mean:
%       p2([envelope(data(1)) slidingam(data(1),10)])
%
%    See also: envelope, slidingavg, slidingrms, slidefun, seizfun

%     Version History:
%        Oct.  5, 2008 - initial version
%        Oct.  7, 2008 - now just an alias to SLIDINGMEAN and SEISFUN
%        Nov. 13, 2008 - update to use SLIDINGAVG and SEIZFUN
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2008 at 21:55 GMT

% todo:

% check nargin
error(nargchk(2,8,nargin))

% check data structure
error(seizchk(data,'dep'))

% alias to other functions
data=seizfun(data,@(x)slidingavg(abs(x),nsamples,varargin{:}));

end
