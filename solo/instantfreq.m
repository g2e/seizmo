function [data]=instantfreq(data)
%INSTANTFREQ    Returns estimated instantaneous frequency of SEIZMO records
%
%    Usage:    data=instantfreq(data)
%
%    Description: INSTANTFREQ(DATA) returns the instantaneous frequency
%     (the point-by-point frequency of a record's analytic signal) of 
%     records in SEIZMO struct DATA.  The instantaneous frequency is
%     estimated from a discrete derivative of the unwrapped instantaneous
%     phase, which is related to the Hilbert transform and the envelope of
%     a signal by the following:
%
%       A(X) = X + iH(X)
%
%                        i*PHI(X)
%       A(X) = ENV(X) * e
%
%     where X is the records' data, i is sqrt(-1), H() denotes a Hilbert
%     transform, A is the analytic signal, ENV is the envelope of the
%     signal and PHI is the instantaneous phase.  The instantaneous
%     frequency is related to the phase by:
%
%     FREQ(X)=d(UNWRAP(PHI(X)))/(2*PI*dX)
%
%     where d stands for differentiation.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Plot up the various components of the analytic signal:
%      plot1([data hilbrt(data) envelope(data) ...
%             instantphase(data) instantfreq(data)])
%
%    See also: HILBRT, ENVELOPE, INSTANTPHASE

%     Version History:
%        June 16, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2010 at 10:50 GMT

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
    [npts,ncmp,delta]=getheader(data,'npts','ncmp','delta');
    nspts=2.^(nextpow2n(npts)+1);
    
    % detail message
    if(verbose)
        disp('Getting Instantaneous Frequency of Record(s)');
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
        
        % get instantaneous frequency (gradient of unwrapped phase)
        % - using discrete approx to derivative
        data(i).dep=qgrad(unwrap(data(i).dep,[],1))/(2*pi*delta(i));
        
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
    error(lasterror)
end

end

function [x]=qgrad(x)
%QGRAD    Quick gradient along 1st dimension

[n,c]=size(x);
if(n>2)
    x=[x(2,:)-x(1,:);
       x(3:n,:)-x(1:n-2,:);
       x(n,:)-x(n-1,:)];
elseif(n>1)
    x=[x(2,:)-x(1,:);
       x(n,:)-x(n-1,:)];
else
    x=zeros([n c],class(x));
end

end

