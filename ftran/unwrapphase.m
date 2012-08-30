function [data]=unwrapphase(data,tol)
%UNWRAPPHASE    Unwraps the phase of SEIZMO records
%
%    Usage:    data=unwrapphase(data)
%              data=unwrapphase(data,tol)
%
%    Description:
%     DATA=UNWRAPPHASE(DATA) unwraps phase data stored in the record(s) of
%     SEIZMO struct DATA by changing absolute jumps (aka phase
%     discontinuities) greater than or equal to PI to their 2*PI
%     complement.  Note that spectral files are converted from Real-
%     Imaginary data format to the Amplitude-Phase data format.  Time
%     Series and XY records may also be "unwrapped" as if their data are
%     phase data.
%
%     DATA=UNWRAPPHASE(DATA,TOL) allows changing the phase jump tolerance
%     TOL from PI (the default) to something else.
%
%    Notes:
%     - Records with IFTYPE 'irlim' are converted to 'iamph'
%     - This is not like SAC's unwrap function, which is basically:
%       data=unwrapphase(dft(data))
%
%    Header changes: IFTYPE, DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Plot the unwrapped phase of time series data:
%     plot2(unwrapphase(keepph(dft(data))));
%
%     % Plot the unwrapped instantaneous phase of time series data:
%     plot2(unwrapphase(instantphase(data)));
%
%    See also: INSTANTPHASE, KEEPPH, GETSPECTRALCMP, DFT

%     Version History:
%        Oct. 20, 2009 - initial version
%        Feb.  3, 2010 - seizmoverbose support, added tol option
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Dec. 21, 2011 - doc update, better checkheader usage
%        Mar. 13, 2012 - use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data,...
        'XYZ_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt unwrapping of phase
try
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % check/expand tol
    if(nargin==1 || isempty(tol)); tol=pi; end
    if(~isreal(tol) || ~any(numel(tol)==[1 nrecs]))
        error('seizmo:unwrapphase:badTol',...
            'TOL must be real scalar or an array w/ 1 real per record!');
    end
    if(isscalar(tol)); tol=tol(ones(nrecs,1),1); end
    
    % get header info
    iftype=getheader(data,'iftype id');
    
    % find spectral
    rlim=strcmpi(iftype,'irlim');
    spectral=rlim | strcmpi(iftype,'iamph');
    
    % detail message
    if(verbose)
        disp('Unwrapping Phase Data of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % convert spectral to amph
    if(any(rlim)); data(rlim)=rlim2amph(data(rlim)); end
    
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
        
        % unwrap
        if(spectral(i))
            % only the phase component of spectral data
            data(i).dep(:,2:2:end)=...
                oclass(unwrap(data(i).dep(:,2:2:end),tol(i),1));
        else
            % time-series and general xy data
            data(i).dep=oclass(unwrap(data(i).dep,tol(i),1));
        end
    
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
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
