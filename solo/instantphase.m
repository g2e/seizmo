function [data]=instantphase(data)
%INSTANTPHASE    Return instantaneous phase of SEIZMO records
%
%    Usage:    data=instantphase(data)
%
%    Description: INSTANTPHASE(DATA) returns the instantaneous phase (the
%     point-by-point phase of a record's analytic signal) of records in
%     SEIZMO struct DATA.  The instantaneous phase is related to the
%     Hilbert transform and the envelope of a signal by the following:
%
%       A(X) = X + iH(X)
%
%                        i*PHI(X)
%       A(X) = ENV(X) * e
%
%     where X is the records' data, i is sqrt(-1), H() denotes a Hilbert
%     transform, A is the analytic signal, ENV is the envelope of the
%     signal and PHI is the instantaneous phase.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Plot up the various components of the analytic signal:
%      plot1([data hilbrt(data) envelope(data) instantphase(data)])
%
%    See also: HILBRT, ENVELOPE, INSTANTFREQ

%     Version History:
%        Oct. 19, 2009 - initial version
%        Jan. 29, 2010 - seizmoverbose support, better warnings
%        June 16, 2010 - code reduction
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2010 at 10:20 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt instant phase
try
    % check header
    data=checkheader(data,'NONTIME_IFTYPE','ERROR','FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % get header info
    [npts,ncmp]=getheader(data,'npts','ncmp');
    nspts=2.^(nextpow2n(npts)+1);
    
    % detail message
    if(verbose)
        disp('Getting Instantaneous Phase of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % frequency domain multiplier to give analytic signal
        h=zeros(nspts(i),1);
        h([1 nspts(i)/2+1])=1;
        h(2:(nspts(i)/2))=2;
        
        % get instantaneous phase
        data(i).dep=angle(ifft(...
            fft(data(i).dep,nspts(i),1).*h(:,ones(1,ncmp(i))),[],1));
        
        % truncate to original length and change class back
        data(i).dep=oclass(data(i).dep(1:npts(i),:));
        
        % dep*
        depmen(i)=mean(data(i).dep(:)); 
        depmin(i)=min(data(i).dep(:)); 
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update header
    data=changeheader(data,...
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
