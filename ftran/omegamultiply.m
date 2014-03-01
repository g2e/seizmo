function [data]=omegamultiply(data)
%OMEGAMULTIPLY    Differentiate SEIZMO records in the frequency domain
%
%    Usage:    data=omegamultiply(data)
%
%    Description:
%     OMEGAMULTIPLY(DATA) basically multiplies each point in the dependent
%     component(s) of spectral records in SEIZMO struct DATA by:
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
%    Header Changes: DEPMEN, DEPMIN, DEPMAX, IDEP
%
%    Examples:
%     % Differentiate spectral data in the time domain vs frequency domain:
%     data=dft(differentiate(idft(data)))
%     data=omegamultiply(data)
%
%    See also: OMEGADIVIDE, OMEGASHIFT, OMEGAANALYTIC, DFT, IDFT,
%              DIFFERENTIATE, INTEGRATE, OMEGAGAUSSIAN

%     Version History:
%        May  12, 2008 - initial version
%        June 11, 2008 - doc cleanup
%        July 19, 2008 - doc update, single ch call, dataless support,
%                        .dep rather than .x
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema (now OMEGAMULTIPLY),
%                        changes idep field
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision),
%                        fixed NPTS bug, full idep support
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        May   6, 2010 - slimmer code for units exchange
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  4, 2012 - doc update, multiplyomega to omegamultiply
%        Feb. 14, 2013 - use strcmpi for consistency
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain differentiation
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
    [delta,npts,idep,iftype]=getheader(data,...
        'delta','npts','idep id','iftype id');
    npts2=npts/2;

    % detail message
    if(verbose)
        disp('Differentiating Record(s) in the Frequency Domain');
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

        % differentiate
        cols=size(data(i).dep,2)/2;
        omega=2*pi*delta(i);
        omega=[0:omega:omega*npts2(i) (npts2(i)-1)*omega:-omega:omega].';
        if(strcmpi(iftype(i),'irlim'))
            % rlim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 0Hz real/imag == 0 else
            % real=-imag*omega & imag=real*omega
            data(i).dep(:,[1:2:end 2:2:end])=oclass(...
                [-data(i).dep(:,2:2:end).*omega(:,ones(1,cols))...
                data(i).dep(:,1:2:end).*omega(:,ones(1,cols))]);
        else % iamph
            % amph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % amp at 0Hz == 0 else amp=amp*omega
            % phase=phase-pi/2 at -Hz and 0Hz else phase=phase+pi/2
            data(i).dep(:,1:2:end)=...
                data(i).dep(:,1:2:end).*omega(:,ones(1,cols));
            data(i).dep(2:(npts2(i)+1),2:2:end)=...
                data(i).dep(2:(npts2(i)+1),2:2:end)+pi/2;
            data(i).dep([1 npts2(i)+2:end],2:2:end)=...
                data(i).dep([1 npts2(i)+2:end],2:2:end)-pi/2;
            data(i).dep=oclass(data(i).dep);
        end

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % change dependent component type
    newunit={'ipop' 'icrackle' 'isnap' 'ijerk' 'iacc' 'ivel' 'idisp' ...
        'iabsmnt' 'iabsity' 'iabseler' 'iabserk' 'iabsnap' 'iabsackle'};
    [tf,idx]=ismember(idep,[newunit(2:end) {'iabspop'}]);
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
    error(lasterror);
end

end
