function [data]=slidingabsmean(data,varargin)
%SLIDINGABSMEAN    Returns sliding-window absolute-mean of SEIZMO records
%
%    Usage:    data=slidingabsmean(data,n)
%              data=slidingabsmean(...,'position','center'|'trail'|'lead')
%              data=slidingabsmean(...,'offset',offset)
%              data=slidingabsmean(...,'edge','truncate'|'pad')
%              data=slidingabsmean(...,'dim',n)
%              data=slidingabsmean(...,'custom',window)
%
%    Description: SLIDINGABSMEAN(DATA,N) applies a centered sliding-window 
%     absolute-mean of 2N+1 samples to the dependent component(s) of 
%     SEIZMO data records in DATA.  N can be a scalar (each record has the
%     same window size) or a vector (define each record's window size 
%     separately).  Sliding windows extending outside the record are
%     truncated (look at 'EDGE' option to change this).
%
%     SLIDINGABSMEAN(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the
%     position of the sliding window relative to the reference data point
%     (the data point which is assigned the window's average).  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1, while for TRAIL or LEAD the window size
%     is N.  Default position is CENTER.
%     
%     SLIDINGABSMEAN(...,'OFFSET',OFFSET) sets the offset of the sliding
%     window from the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of 
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     OFFSET may be a vector of offsets specifying each record's offset.
%     
%     SLIDINGABSMEAN(...,'EDGE','TRUNCATE'|'PAD') sets the handling of edge
%     cases.  TRUNCATE eliminates points in the sliding-window that do not 
%     reference a datapoint (for instance if the window extends before or 
%     after the data, that portion of the window will be truncated).  PAD 
%     adds zeros to the data so that all the points of the sliding-window 
%     always reference some value.  This will create a tapered look at the
%     edges of the data.  The default setting of EDGE is TRUNCATE.
%
%     SLIDINGABSMEAN(...,'DIM',N) slides across dimension N raher than the
%     default 1 (1 slides down the component - 2 would slide across the
%     components).
%
%     SLIDINGABSMEAN(...,'CUSTOM',WINDOW) allows a custom sliding window
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
%     - SLIDINGABSMEAN is faster than SLIDEFUN because it uses SLIDINGAVG
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%      Compare an envelope and a 21-sample sliding-window absolute-mean:
%       p2([envelope(data(1)) slidingabsmean(data(1),10)])
%
%    See also: ENVELOPE, SLIDINGAVG, SLIDINGRMS, SLIDINGFUN, SEIZMOFUN

%     Version History:
%        Oct.  5, 2008 - initial version
%        Oct.  7, 2008 - now just an alias to SLIDINGMEAN and SEISFUN
%        Nov. 13, 2008 - update to use SLIDINGAVG and SEIZFUN
%        Nov. 22, 2008 - update for new name schema (now SLIDINGABSMEAN)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June  9, 2009 - nsamples now in varargin, toggle seizmocheck,
%                        up nargin allowed
%        Oct. 13, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 13, 2009 at 06:05 GMT

% todo:

% check nargin
msg=nargchk(2,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% alias to other functions
data=seizmofun(data,@(x)slidingavg(abs(x),varargin{:}));

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
