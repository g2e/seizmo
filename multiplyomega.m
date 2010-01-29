function [data]=multiplyomega(data)
%MULTIPLYOMEGA    Differentiate SEIZMO records in the frequency domain
%
%    Usage:    data=multiplyomega(data)
%
%    Description: MULTIPLYOMEGA(DATA) basically multiplies each point in
%     the dependent component(s) of spectral files by:
%       OMEGA=2.0 * PI * FREQ
%     to perform the equivalent of differentiation in the time domain.  
%     This is particularly handy when working with spectral data as it 
%     avoids the forward and inverse fourier transform necessary for time 
%     domain differentiation.
%
%    Notes:
%     - Read the source code below for a better description of the
%       operations performed for frequency-domain differentiation.
%
%    Header Changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Differentiate spectral data in the time domain vs frequency domain:
%      data=dft(differentiate(idft(data)))
%      data=multiplyomega(data)
%
%    See also: DIVIDEOMEGA, DFT, IDFT

%     Version History:
%        May  12, 2008 - initial version
%        June 11, 2008 - doc cleanup
%        July 19, 2008 - doc update, single ch call, dataless support,
%                        .dep rather than .x
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema (now MULTIPLYOMEGA),
%                        changes idep field
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision),
%                        fixed NPTS bug, full idep support
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2010 at 18:55 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain differentiation
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
        error('seizmo:multiplyomega:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    elseif(any(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph')))
        error('seizmo:multiplyomega:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph'))) ...
            '\nDatatype of record(s) in DATA must be spectral!']);
    end

    % detail message
    if(verbose)
        disp('Differentiating Record(s) in the Frequency Domain');
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

        % differentiate
        cols=size(data(i).dep,2);
        omega=[0 2*pi*[(1:npts2(i))*delta(i) ...
            ((npts2(i)-1):-1:1)*delta(i)]].';
        if(strcmp(iftype(i),'irlim'))
            % rlim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 0Hz real/imag == 0 else
            % real=-imag*omega & imag=real*omega
            data(i).dep(:,[1:2:end 2:2:end])=oclass(...
                [-data(i).dep(:,2:2:end).*omega(:,ones(1,cols/2))...
                data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2))]);
        else % iamph
            % amph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % amp at 0Hz == 0 else amp=amp*omega
            % phase=phase-pi/2 at -Hz and 0Hz else phase=phase+pi/2
            data(i).dep(:,1:2:end)=...
                data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2));
            data(i).dep(2:(npts2(i)+1),2:2:end)=...
                data(i).dep(2:(npts2(i)+1),2:2:end)+pi/2;
            data(i).dep([1 npts2(i)+2:end],2:2:end)=...
                data(i).dep([1 npts2(i)+2:end],2:2:end)-pi/2;
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
    iscrackle=strcmpi(idep,'icrackle');
    issnap=strcmpi(idep,'isnap');
    isjerk=strcmpi(idep,'ijerk');
    isacc=strcmpi(idep,'iacc');
    isvel=strcmpi(idep,'ivel');
    isdisp=strcmpi(idep,'idisp');
    isabsmnt=strcmpi(idep,'iabsmnt');
    isabsity=strcmpi(idep,'iabsity');
    isabseler=strcmpi(idep,'iabseler');
    isabserk=strcmpi(idep,'iabserk');
    isabsnap=strcmpi(idep,'iabsnap');
    isabsackl=strcmpi(idep,'iabsackl');
    isabspop=strcmpi(idep,'iabspop');
    idep(iscrackle)={'ipop'};
    idep(issnap)={'icrackle'};
    idep(isjerk)={'isnap'};
    idep(isacc)={'ijerk'};
    idep(isvel)={'iacc'};
    idep(isdisp)={'ivel'};
    idep(isabsmnt)={'idisp'};
    idep(isabsity)={'iabsmnt'};
    idep(isabseler)={'iabsity'};
    idep(isabserk)={'iabseler'};
    idep(isabsnap)={'iabserk'};
    idep(isabsackl)={'iabsnap'};
    idep(isabspop)={'iabsackl'};
    idep(~(iscrackle | issnap | isjerk | isacc | isvel | isdisp | isabsmnt |...
        isabsity | isabseler | isabserk | isabsnap | isabsackl | isabspop))...
        ={'iunkn'};

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
