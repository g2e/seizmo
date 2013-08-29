function [data]=slidingmean(data,n,varargin)
%SLIDINGMEAN    Returns sliding-window mean of SEIZMO records
%
%    Usage:    data=slidingmean(data,n)
%              data=slidingmean(...,'position','center'|'trail'|'lead')
%              data=slidingmean(...,'offset',offset)
%              data=slidingmean(...,'edge','truncate'|'pad')
%              data=slidingmean(...,'nans','skipped'|'too'|'only')
%              data=slidingmean(...,'dim',n)
%              data=slidingmean(...,'custom',window)
%
%    Description:
%     DATA=SLIDINGMEAN(DATA,N) applies a centered sliding-window mean of
%     2N+1 samples to the dependent component(s) of records in SEIZMO
%     struct DATA.  N can be a scalar (each record has the same window
%     size) or a vector (define each record's window size separately).
%     Sliding windows extending outside the record are truncated (look at
%     'EDGE' option to change this).
%
%     DATA=SLIDINGMEAN(...,'POSITION','CENTER'|'TRAIL'|'LEAD')
%     DATA=SLIDINGMEAN(...,'OFFSET',OFFSET)
%     DATA=SLIDINGMEAN(...,'EDGE','TRUNCATE'|'PAD')
%     DATA=SLIDINGMEAN(...,'NANS','SKIPPED'|'TOO'|'ONLY')
%     DATA=SLIDINGMEAN(...,'DIM',N)
%     DATA=SLIDINGMEAN(...,'CUSTOM',WINDOW)
%      See SLIDINGAVG for details!
%
%    Notes:
%     - SLIDINGMEAN is faster than SLIDEFUN because it uses SLIDINGAVG
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Compare a 21-sample sliding-window mean to the original record:
%     plot2([data(1) slidingmean(data(1),10)])
%
%    See also: SLIDINGABSMEAN, SLIDINGAVG, SLIDINGRMS, SOLOFUN

%     Version History:
%        Feb.  3, 2010 - initial version
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

% attempt sliding mean
% NOTE: b/c slidingmean allows setting n per record we cannot reduce this
%       function to the single solofun call:
%        data=solofun(data,@(x)slidingavg(x,n,varargin{:}));
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check n
    if(~isreal(n) || ~any(numel(n)==[1 nrecs]) || any(n<0))
        error('seizmo:slidingmean:badInput',...
            'N must be a positive real scalar or vector!');
    end
    if(isscalar(n)); n(1:nrecs,1)=n; end
    
    % detail message
    if(verbose)
        disp('Appyling Sliding Mean to Dependent Data of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % apply function to records
    ncmp=nan(nrecs,1); npts=ncmp;
    depmen=ncmp; depmin=ncmp; depmax=ncmp;
    for i=1:nrecs
        oclass=str2func(class(data(i).dep));
        data(i).dep=oclass(...
            slidingavg(double(data(i).dep),n(i),varargin{:}));
        
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
