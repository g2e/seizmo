function [data]=omegaanalytic(data)
%OMEGAANALYTIC    Analytic signal of SEIZMO records in the frequency domain
%
%    Usage:    data=omegaanalytic(data)
%
%    Description:
%     DATA=OMEGAANALYTIC(DATA) returns the analytic signal of SEIZMO
%     spectral records in DATA.  The analytic signal is useful for
%     determining the envelope and instantaneous phase of a signal.
%
%    Notes:
%     - IDFT needs to be passed the 'nonsymmetric' flag in order to
%       preserve the complex values of the analytic signal when converting
%       to the time domain!
%     - Because DFT zero-pads records this operation does not give the same
%       result as Matlab's hilbert.
%
%    Examples:
%     % Get the envelope and instantaneous phase:
%     data=idft(omegaanalytic(data),'nonsymmetric');
%     edata=solofun(data,@abs);
%     pdata=solofun(data,@angle);
%
%    See also: OMEGAHILBERT, OMEGASHIFT, OMEGADIVIDE, OMEGAMULTIPLY,
%              OMEGAGAUSSIAN, DFT, IDFT, HILBRT, ENVELOPE, INSTANTPHASE,
%              INSTANTFREQ

%     Version History:
%        Feb.  5, 2012 - initial version
%        Feb. 14, 2013 - use strcmpi for consistency
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain analytic signal
try
    % check headers
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % retreive header info
    [npts,iftype]=getheader(data,'npts','iftype id');
    npts2=npts/2;
    
    % detail message
    if(verbose)
        disp('Getting Record(s) Analytic Signal in the Frequency Domain');
        print_time_left(0,nrecs);
    end

    % loop through records
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

        % analytic signal
        if(strcmpi(iftype(i),'irlim'))
            data(i).dep(2:npts2(i),:)=data(i).dep(2:npts2(i),:)*2;
            data(i).dep(npts2(i)+2:end,:)=0;
        else % iamph
            data(i).dep(2:npts2(i),1:2:end)=...
                data(i).dep(2:npts2(i),1:2:end)*2;
            data(i).dep(npts2(i)+2:end,1:2:end)=0;
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
