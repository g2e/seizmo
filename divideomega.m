function [data]=divideomega(data)
%DIVIDEOMEGA    Integrate SEIZMO records in the frequency domain
%
%    Usage:    data=divideomega(data)
%
%    Description: DIVIDEOMEGA(DATA) basically divides each point in the 
%     dependent component(s) of spectral files by:
%       OMEGA=2.0 * PI * FREQ
%     to perform the equivalent of integration in the time domain.  This is
%     particularly handy when working with spectral data as it avoids
%     the forward and inverse fourier transform necessary for time domain 
%     integration.  It is also useful for reducing the dynamic range of 
%     spectral data.
%
%    Notes:
%     - Read the source code below for a better description of the
%       operations performed for frequency-domain integration.
%
%    Header Changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Integrate spectral data in the time domain vs frequency domain:
%      data=dft(integrate(idft(data)))
%      data=divideomega(data)
%
%    See also: MULTIPLYOMEGA, DFT, IDFT

%     Version History:
%        May  12, 2008 - initial version
%        June 11, 2008 - doc cleanup
%        July  8, 2008 - doc update, single ch call, .dep rather than .x
%        July 19, 2008 - doc update, dataless support
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema (now DIVIDEOMEGA),
%                        changes idep field
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision),
%                        fixed NPTS bug, full idep support
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        May   6, 2010 - slimmer code for units exchange
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   6, 2010 at 23:10 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain integration
try
    % check headers
    data=checkheader(data);

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % retreive header info
    leven=getlgc(data,'leven');
    [idep,iftype]=getenumid(data,'idep','iftype');
    [e,delta,npts]=getheader(data,'e','delta','npts');
    npts2=npts/2;

    % check leven,iftype
    if(any(strcmpi(leven,'false')))
        error('seizmo:divideomega:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    elseif(any(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph')))
        error('seizmo:divideomega:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph'))) ...
            '\nDatatype of record(s) in DATA must be spectral!']);
    end
    
    % detail message
    if(verbose)
        disp('Integrating Record(s) in the Frequency Domain');
        print_time_left(0,nrecs);
    end

    % loop through records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
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

        % integrate
        cols=size(data(i).dep,2);
        omega=[0 1./(2*pi*[(1:npts2(i))*delta(i) ...
            ((npts2(i)-1):-1:1)*delta(i)])].';
        if(strcmp(iftype(i),'irlim'))
            % rlim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 0Hz real/imag == 0 else
            % real=imag/omega & imag=-real/omega
            data(i).dep(:,[1:2:end 2:2:end])=oclass(...
                [data(i).dep(:,2:2:end).*omega(:,ones(1,cols/2)) ...
                -data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2))]);
        else % iamph
            % amph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % amp at 0Hz == 0 else amp=amp/omega
            % phase=phase+pi/2 at -Hz and 0Hz else phase=phase-pi/2
            data(i).dep(:,1:2:end)=...
                data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2));
            data(i).dep(2:(npts2(i)+1),2:2:end)=...
                data(i).dep(2:(npts2(i)+1),2:2:end)-pi/2;
            data(i).dep([1 npts2(i)+2:end],2:2:end)=...
                data(i).dep([1 npts2(i)+2:end],2:2:end)+pi/2;
            data(i).dep=oclass(data(i).dep);
        end

        % dep*
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % change dependent component type
    newunit={'icrackle' 'isnap' 'ijerk' 'iacc' 'ivel' 'idisp' 'iabsmnt' ...
        'iabsity' 'iabseler' 'iabserk' 'iabsnap' 'iabsackle' 'iabspop'};
    [tf,idx]=ismember(idep,[{'ipop'} newunit(1:end-1)]);
    idep(tf)=newunit(idx(tf));
    idep(~tf)={'iunkn'};

    % update header
    data=changeheader(data,'idep',idep,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
