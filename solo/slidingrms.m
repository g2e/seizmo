function [data]=slidingrms(data,n,varargin)
%SLIDINGRMS    Returns sliding-window root-mean-square of SEIZMO records
%
%    Usage:    data=slidingrms(data,n)
%              data=slidingrms(...,'position','center'|'trail'|'lead')
%              data=slidingrms(...,'offset',offset)
%              data=slidingrms(...,'edge','truncate'|'pad')
%              data=slidingrms(...,'dim',n)
%              data=slidingrms(...,'custom',window)
%
%    Description:
%     DATA=SLIDINGRMS(DATA,N) applies a centered sliding-window rms
%     (root-mean-square) of 2N+1 samples to the dependent component(s) of 
%     SEIZMO data records in DATA.  N can be a scalar (each record has the
%     same window size) or a vector (define each record's window size 
%     separately).  Sliding windows extending outside the record are
%     truncated (look at 'EDGE' option to change this).
%
%     DATA=SLIDINGRMS(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the
%     position of the sliding window relative to the reference data point.
%     CENTER positions the window such that the reference point is at its
%     center.  TRAIL positions the window to trail the reference point such
%     that the reference point has the highest index in the window.  LEAD
%     sets the window to lead the reference point such that the reference
%     point has the lowest index in the window.  Note that the window size
%     for a CENTER positioning is 2N+1, while for TRAIL or LEAD the window
%     size is N.  Default position is CENTER.
%     
%     DATA=SLIDINGRMS(...,'OFFSET',OFFSET) sets the offset of the sliding
%     window from to the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of 
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     OFFSET may be a vector of offsets specifying each record's offset.
%     
%     DATA=SLIDINGRMS(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the
%     edge cases.  TRUNCATE eliminates points in the sliding-window that do
%     not reference a data point (ie if the window extends before or after
%     the data, that portion of the window will be truncated).  PAD adds
%     zeros to the data so that all the points of the sliding-window always
%     reference some value.  Default setting is TRUNCATE.
%
%     DATA=SLIDINGRMS(...,'DIM',N) specifies an alternative dimension to
%     slide across rather than the default of 1 (N=1 slides down the
%     component, N=2 slides across the components).
%
%     DATA=SLIDINGRMS(...,'CUSTOM',WINDOW) allows a custom sliding window
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
%     - SLIDINGRMS is faster than SLIDINGFUN because it uses SLIDINGAVG
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Compare an envelope and 21-sample sliding-window root-mean-square:
%     plot2([envelope(data(1)) slidingrms(data(1),10)])
%
%    See also: ENVELOPE, SLIDINGABSMEAN, SLIDINGFUN, SLIDINGAVG, SOLOFUN

%     Version History:
%        Apr.  9, 2008 - initial version
%        Apr. 23, 2008 - changed behavior for windows with even npts
%        May  12, 2998 - dep* fix
%        July 17, 2008 - history update, doc update, .dep rather
%                        than .x, dataless handling
%        Oct.  5, 2008 - big change: updated to match options of SLIDEFUN,
%                        name changed from RMS to SLIDINGRMS
%        Oct.  7, 2008 - now just an alias to SLIDINGMEAN and SEISFUN
%        Nov. 13, 2008 - update to use SLIDINGAVG and SEIZFUN
%        Nov. 22, 2008 - update for new names again
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June  9, 2009 - nsamples now in varargin, toggle seizmocheck,
%                        up nargin allowed
%        Feb.  3, 2010 - proper SEIZMO handling, versioninfo caching
%        Apr. 22, 2010 - allow multiple N as advertised
%        Jan.  6, 2011 - drop versioninfo caching, nargchk fix,
%                        seizmofun/solofun rename
%        Apr.  3, 2012 - minor doc update
%        May  30, 2012 - allow N=0
%        May  31, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2012 at 09:45 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt sliding rms
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check n
    if(~isreal(n) || ~any(numel(n)==[1 nrecs]) || any(n<0))
        error('seizmo:slidingrms:badInput',...
            'N must be a positive real scalar or vector!');
    end
    if(isscalar(n)); n(1:nrecs,1)=n; end
    
    % detail message
    if(verbose)
        disp('Appyling Sliding RMS to Dependent Data of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % apply function to records
    ncmp=nan(nrecs,1); npts=ncmp;
    depmen=ncmp; depmin=ncmp; depmax=ncmp;
    for i=1:nrecs
        oclass=str2func(class(data(i).dep));
        data(i).dep=oclass(...
            sqrt(slidingavg(double(data(i).dep).^2,n(i),varargin{:})));
        
        % get npts, ncmp, dep*
        [npts(i),ncmp(i)]=size(data(i).dep);
        if(npts(i)) % skip dep* for empty
            depmen(i)=nanmean(data(i).dep(:)); 
            depmin(i)=min(data(i).dep(:)); 
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update header
    data=changeheader(data,'npts',npts,'ncmp',ncmp,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
