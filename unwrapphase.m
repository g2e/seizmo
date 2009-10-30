function [data]=unwrapphase(data)
%UNWRAPPHASE    Unwraps the phase of SEIZMO records
%
%    Usage:    data=unwrapphase(data)
%
%    Description: DATA=UNWRAPPHASE(DATA) unwraps phase data stored in
%    SEIZMO struct DATA by changing absolute jumps greater than or equal to
%    PI to their 2*PI complement.  Note that spectral files are converted
%    from Real-Imaginary data format to the Amplitude-Phase data format in
%    this process.
%
%    Notes:
%     - Records with IFTYPE 'irlim' are converted to 'iamph'
%     - This is not like SAC's unwrap function, which is basically:
%       data=unwrapphase(dft(data))
%
%    Header changes: IFTYPE, DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Plot the unwrapped phase of time series data:
%      plot1(unwrapphase(keepph(dft(data)));
%
%     Plot the unwrapped instantaneous phase of time series data:
%      plot1(unwrapphase(instantphase(data)));
%
%    See also: INSTANTPHASE, KEEPPH, GETSPECTRALCMP, DFT

%     Version History:
%        Oct. 20, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 20, 2009 at 23:45 GMT

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
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt unwrapping of phase
try
    % get header info
    iftype=getenumid(data,'iftype');
    
    % find spectral
    rlim=strcmpi(iftype,'irlim');
    spectral=rlim | strcmpi(iftype,'iamph');
    
    % cannot do xyz records
    if(any(strcmpi(iftype,'ixyz')))
        error('seizmo:instantphase:badIFTYPE',...
            'Illegal operation on XYZ data!');
    end
    
    % convert spectral to amph
    if(any(rlim)); data(rlim)=rlim2amph(data(rlim)); end
    
    % loop over records
    nrecs=numel(data);
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep)); continue; end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % unwrap
        if(spectral(i))
            % only the phase component of spectral data
            data(i).dep(:,2:2:end)=...
                oclass(unwrap(data(i).dep(:,2:2:end),[],1));
        else
            % time-series and general xy data
            data(i).dep=oclass(unwrap(data(i).dep,[],1));
        end
    
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
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
