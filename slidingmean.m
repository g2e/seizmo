function [data]=slidingmean(data,varargin)
%SLIDINGMEAN    Returns sliding-window mean of SEIZMO records
%
%    Usage:    data=slidingmean(data,n)
%              data=slidingmean(...,'position','center'|'trail'|'lead')
%              data=slidingmean(...,'offset',offset)
%              data=slidingmean(...,'edge','truncate'|'pad')
%              data=slidingmean(...,'dim',n)
%              data=slidingmean(...,'custom',window)
%
%    Description: SLIDINGMEAN(DATA,N) applies a centered sliding-window 
%     mean of 2N+1 samples to the dependent component(s) of records in
%     SEIZMO struct DATA.  N can be a scalar (each record has the same
%     window size) or a vector (define each record's window size 
%     separately).  Sliding windows extending outside the record are
%     truncated (look at 'EDGE' option to change this).
%
%     SLIDINGMEAN(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the
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
%     SLIDINGMEAN(...,'OFFSET',OFFSET) sets the offset of the sliding
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
%     SLIDINGMEAN(...,'EDGE','TRUNCATE'|'PAD') sets the handling of edge
%     cases.  TRUNCATE eliminates points in the sliding-window that do not 
%     reference a datapoint (for instance if the window extends before or 
%     after the data, that portion of the window will be truncated).  PAD 
%     adds zeros to the data so that all the points of the sliding-window 
%     always reference some value.  This will create a tapered look at the
%     edges of the data.  The default setting of EDGE is TRUNCATE.
%
%     SLIDINGMEAN(...,'DIM',N) slides across dimension N raher than the
%     default 1 (1 slides down the component - 2 would slide across the
%     components).
%
%     SLIDINGMEAN(...,'CUSTOM',WINDOW) allows a custom sliding window
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
%     - SLIDINGMEAN is faster than SLIDEFUN because it uses SLIDINGAVG
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%      Compare a 21-sample sliding-window mean to the original record:
%       p2([data(1) slidingmean(data(1),10)])
%
%    See also: SLIDINGABSMEAN, SLIDINGAVG, SLIDINGRMS, SEIZMOFUN

%     Version History:
%        Feb.  3, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2010 at 17:05 GMT

% todo:

% check nargin
msg=nargchk(2,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt sliding rms
try
    % alias to other functions
    data=seizmofun(data,@(x)slidingavg(x,varargin{:}));
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
