function [data]=hilbrt(data)
%HILBRT    Return Hilbert transform of SEIZMO records
%
%    Usage:    data=hilbrt(data)
%
%    Description:
%     HILBRT(DATA) calculates and returns the Hilbert transform of SEIZMO
%     records.  The Hilbert tranform is basically a -pi/2 phase shift.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Quick way to do a positive pi/2 phase shift:
%     data=hilbrt(multiply(data,-1))
%
%    See also: ENVELOPE, INSTANTPHASE, INSTANTFREQ, OMEGAHILBERT,
%              OMEGAANALYTIC

%     Version History:
%        Jan. 30, 2008 - initial version
%        Feb. 28, 2008 - seischk support and class support
%        Mar.  4, 2008 - doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - minor doc update
%        July 17, 2008 - history update, doc update, now uses
%                        SEIZMO functions rather than Matlab's hilbert
%        Nov. 22, 2008 - doc update, update for new name schema
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Oct. 15, 2009 - does data and header checking now, try/catch will
%                        assure checking state is kept correct, no longer
%                        calls other SEIZMO functions
%        Oct. 19, 2009 - added checks for IFTYPE and LEVEN
%        Jan. 29, 2010 - seizmoverbose support, better warnings
%        June 16, 2010 - code reduction
%        Feb.  4, 2012 - no zero-padding, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  4, 2012 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt hilbert
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
    [npts,ncmp]=getheader(data,'npts','ncmp');
    
    % detail message
    if(verbose)
        disp('Getting Hilbert Transform of Record(s)');
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
        
        % get hilbert
        data(i).dep=oclass(imag(ifft(...
            fft(data(i).dep,[],1).*h(:,ones(1,ncmp(i))),[],1)));
        
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
