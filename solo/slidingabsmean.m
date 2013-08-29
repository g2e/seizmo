function [data]=slidingabsmean(data,n,varargin)
%SLIDINGABSMEAN    Returns sliding-window absolute-mean of SEIZMO records
%
%    Usage:    data=slidingabsmean(data,n)
%              data=slidingabsmean(...,'position','center'|'trail'|'lead')
%              data=slidingabsmean(...,'offset',offset)
%              data=slidingabsmean(...,'edge','truncate'|'pad')
%              data=slidingabsmean(...,'nans','skipped'|'too'|'only')
%              data=slidingabsmean(...,'dim',n)
%              data=slidingabsmean(...,'custom',window)
%
%    Description:
%     DATA=SLIDINGABSMEAN(DATA,N) applies a centered sliding-window 
%     absolute-mean of 2N+1 samples to the dependent component(s) of 
%     SEIZMO data records in DATA.  N can be a scalar (each record has the
%     same window size) or a vector (define each record's window size 
%     separately).  Sliding windows extending outside the record are
%     truncated (look at 'EDGE' option to change this).
%
%     DATA=SLIDINGABSMEAN(...,'POSITION','CENTER'|'TRAIL'|'LEAD')
%     DATA=SLIDINGABSMEAN(...,'OFFSET',OFFSET)
%     DATA=SLIDINGABSMEAN(...,'EDGE','TRUNCATE'|'PAD')
%     DATA=SLIDINGABSMEAN(...,'NANS','SKIPPED'|'TOO'|'ONLY')
%     DATA=SLIDINGABSMEAN(...,'DIM',N)
%     DATA=SLIDINGABSMEAN(...,'CUSTOM',WINDOW)
%      See SLIDINGAVG for details!
%
%    Notes:
%     - SLIDINGABSMEAN is faster than SLIDEFUN because it uses SLIDINGAVG
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Compare an envelope and a 21-sample sliding-window absolute-mean:
%     plot2([envelope(data(1)) slidingabsmean(data(1),10)])
%
%    See also: ENVELOPE, SLIDINGAVG, SLIDINGRMS, SLIDINGFUN, SOLOFUN

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
%        Feb.  3, 2010 - proper SEIZMO handling, versioninfo caching
%        Apr. 22, 2010 - allow multiple N as advertised
%        Jan.  6, 2011 - drop versioninfo caching, nargchk fix,
%                        seizmofun/solofun rename
%        Apr.  3, 2012 - minor doc update
%        May  30, 2012 - allow N=0
%        May  31, 2012 - minor doc update
%        Aug.  2, 2013 - removed most doc repetition for easy maintenance
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  2, 2013 at 09:45 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt sliding absolute mean
% NOTE: b/c slidingabsmean allows setting n per record we cannot reduce
%       this function to the single solofun call:
%        data=solofun(data,@(x)slidingavg(abs(x),n,varargin{:}));
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check n
    if(~isreal(n) || ~any(numel(n)==[1 nrecs]) || any(n<0))
        error('seizmo:slidingabsmean:badInput',...
            'N must be a positive real scalar or vector!');
    end
    if(isscalar(n)); n(1:nrecs,1)=n; end
    
    % detail message
    if(verbose)
        disp('Appyling Sliding Abs. Mean to Dependent Data of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % apply function to records
    ncmp=nan(nrecs,1); npts=ncmp;
    depmen=ncmp; depmin=ncmp; depmax=ncmp;
    for i=1:nrecs
        oclass=str2func(class(data(i).dep));
        data(i).dep=oclass(...
            slidingavg(abs(double(data(i).dep)),n(i),varargin{:}));
        
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
