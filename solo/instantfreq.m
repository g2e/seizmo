function [data]=instantfreq(data)
%INSTANTFREQ    Returns instantaneous frequency of SEIZMO records
%
%    Usage:    data=instantfreq(data)
%
%    Description:
%     INSTANTFREQ(DATA) returns the instantaneous frequency (the point-by-
%     point frequency of a record's analytic signal) of records in SEIZMO
%     struct DATA.  The instantaneous frequency is estimated from a
%     discrete derivative of the unwrapped instantaneous phase, which is
%     related to the Hilbert transform, the analytic signal and the
%     envelope of a signal by the following:
%
%      (1) A(t) = X(t) + iH(X(t))
%
%                           i*PHI(t)
%      (2) A(t) = ENV(t) * e
%
%     where X is the records' data, i is sqrt(-1), H() denotes a Hilbert
%     transform, A is the analytic signal, ENV is the envelope of the
%     signal and PHI is the instantaneous phase.  The instantaneous
%     frequency is related to the instantaneous phase by:
%
%      (3) FREQ(t)=d(UNWRAP(PHI(t)))/(2*PI*dt)
%
%     where d stands for discrete differentiation.
%
%    Notes:
%     - Make sure to remove the trend/mean beforehand!
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Plot up the various components of the analytic signal:
%     plot1([data hilbrt(data) envelope(data) ...
%         instantphase(data) instantfreq(data)])
%
%    See also: HILBRT, ENVELOPE, INSTANTPHASE, OMEGAHILBERT,
%              OMEGAANALYTIC

%     Version History:
%        June 16, 2010 - initial version
%        Feb.  4, 2012 - doc update, no zero-padding
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2010 at 10:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt instant phase
try
    % check header
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % get header info
    [npts,ncmp,delta]=getheader(data,'npts','ncmp','delta');
    
    % detail message
    if(verbose)
        disp('Getting Instantaneous Frequency of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    [depmin,depmen,depmax]=deal(nan(nrecs,1));
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
        % = (1+sgn(w))
        h=zeros(npts(i),1);
        if(mod(npts,2)) % odd
            h(1)=1;
            h(2:((npts(i)+1)/2))=2;
        else % even
            h([1 npts(i)/2+1])=1;
            h(2:(npts(i)/2))=2;
        end
        
        % get instantaneous phase
        data(i).dep=oclass(angle(ifft(...
            fft(data(i).dep,[],1).*h(:,ones(1,ncmp(i))),[],1)));
        
        % get instantaneous frequency (gradient of unwrapped phase)
        % - using discrete approximation to derivative
        data(i).dep=qgrad(unwrap(data(i).dep,[],1))/(2*pi*delta(i));
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:)); 
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

