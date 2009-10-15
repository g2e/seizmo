function [data]=envelope(data)
%ENVELOPE    Return envelopes of SEIZMO records
%
%    Usage:    data=envelope(data)
%
%    Description: ENVELOPE(DATA) returns the envelope (the complex
%     magnitude of a record's analytic signal) of records in SEIZMO struct
%     DATA.  The envelope formula is as follows:
%
%       ENVELOPE = sqrt(X^2 + H(X)^2)
%
%     where X is the records' data and H denotes a Hilbert transform.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Plot the envelopes against the data:
%      recordsection([data; envelope(data)])
%
%    See also: HILBRT

%     Version History:
%        Jan. 30, 2008 - initial version
%        Feb. 23, 2008 - seischk support and class support
%        Feb. 28, 2008 - support for changing fft zeropadding
%        Mar.  4, 2008 - doc update
%        May  12, 2998 - dep* fix
%        July 17, 2008 - history update, doc update, now uses
%                        SEIZMO functions rather than Matlab's hilbert
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Oct. 15, 2009 - does data and header checking now, try/catch will
%                        assure checking state is kept correct, no longer
%                        calls other SEIZMO functions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2009 at 16:20 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    rethrow(lasterror)
end

% attempt envelope
try
    % get header info
    [npts,ncmp]=getheader(data,'npts','ncmp');
    nspts=2.^(nextpow2n(npts)+1);
    
    % loop over records
    nrecs=numel(data);
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep)); continue; end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % multiplier to give analytic signal
        h=zeros(nspts(i),1);
        h([1 nspts(i)/2+1])=1;
        h(2:(nspts(i)/2))=2;
        
        % get envelope
        data(i).dep=abs(ifft(...
            fft(data(i).dep,nspts(i),1).*h(:,ones(1,ncmp(i))),...
            [],1));
        
        % truncate to original length and change class back
        data(i).dep=oclass(data(i).dep(1:npts(i),:));
        
        % dep*
        depmen(i)=mean(data(i).dep(:)); 
        depmin(i)=min(data(i).dep(:)); 
        depmax(i)=max(data(i).dep(:));
    end
    
    % update header
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    rethrow(lasterror)
end

end
