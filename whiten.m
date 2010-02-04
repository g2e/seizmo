function [data]=whiten(data,halfwindow,varargin)
%WHITEN    Spectral whitening/normalization of SEIZMO data records
%
%    Usage:    data=whiten(data)
%              data=whiten(data,halfwindow)
%              data=whiten(...,'optionname',optionvalue,...)
%
%    Description: WHITEN(DATA) will perform spectral normalization (aka
%     whitening) on records in the SEIZMO structure DATA.  Normalization is
%     performed by dividing the complex spectrum by a smoothed version of
%     the amplitude spectrum.  Smoothing utilizes a 41-sample sliding mean.
%     This is NOT equivalent to the 'whiten' command in SAC (see function
%     PREWHITEN).  This operation is particularly suited for ambient noise
%     studies.
%
%     WHITEN(DATA,N) allows changing the half-width of the smoothing window
%     (the width is actually 2N+1).  The default N is 20.  See the function
%     SLIDINGMEAN for more information.
%
%     WHITEN(...,'OPTIONNAME',OPTIONVALUE,...) will pass sliding average
%     options on to the SLIDINGMEAN call.  See SLIDINGMEAN for more
%     information.
%
%    Notes:
%     - Suggested Reading:
%       - Bensen et al, 2007, Processing Seismic Ambient Noise Data to
%         Obtain Reliable Broad-Band Surface Wave Dispersion Measurements,
%         GJI, Vol. 169, p. 1239-1260
%
%    Header Changes: DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     Spectral normalization returns much whiter noise:
%      plot1([data(1) whiten(data(1))])
%
%    See also: SLIDINGMEAN, PREWHITEN, UNPREWHITEN

%     Version History:
%        June  9, 2009 - initial version
%        June 11, 2009 - updated default halfwindow from 2 to 20
%        June 24, 2009 - now transparent to filetype (except ixyz)
%        Dec.  4, 2009 - fixed rlim handling
%        Dec.  7, 2009 - no divide by zero by adding eps to smooth spectra
%        Feb.  3, 2010 - proper SEIZMO handling
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2010 at 18:25 GMT

% todo:

% check number of inputs
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt spectral whitening
try
    % verbosity
    verbose=seizmoverbose;
    
    % detail message
    if(verbose); disp('Beginning Spectral Whitening of Record(s)'); end
    
    % get some header fields
    leven=getlgc(data,'leven');
    iftype=getenumid(data,'iftype');

    % cannot do xyz records
    if(any(strcmpi(iftype,'ixyz')))
        error('seizmo:whiten:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ', ...
            find(~strcmpi(iftype,'ixyz'))) ...
            '\nInvalid operation on XYZ record(s)!']);
    end

    % cannot do unevenly sampled records
    if(any(strcmpi(leven,'false')))
        error('seizmo:whiten:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    end

    % default centered window half size
    if(nargin==1 || isempty(halfwindow)); halfwindow=20; end

    % get filetype logical arrays
    istime=strcmpi(iftype,'itime');
    isxy=strcmpi(iftype,'ixy');
    isrlim=strcmpi(iftype,'irlim');
    isamph=strcmpi(iftype,'iamph');

    % get amph and rlim type records
    if(any(istime | isxy))
        amph(istime | isxy)=dft(data(istime | isxy));
        data(istime | isxy)=amph2rlim(amph(istime | isxy));
    end
    if(any(isrlim))
        amph(isrlim)=rlim2amph(data(isrlim));
    end
    if(any(isamph))
        amph(isamph)=data(isamph);
        data(isamph)=amph2rlim(data(isamph));
    end

    % fake amph records as rlim (to get by dividerecords checks/fixes)
    amph=changeheader(amph,'iftype','irlim');

    % get smoothed amplitude records
    amph=slidingmean(amph,halfwindow,varargin{:});

    % copy amplitude over phase & add eps to avoid divide by zero
    amph=seizmofun(amph,@(x)x(:,[1:2:end; 1:2:end])+eps);

    % divide complex by smoothed amplitude
    data=dividerecords(data,amph);

    % convert back to original type
    if(any(isamph))
        data(isamph)=rlim2amph(data(isamph));
    end
    if(any(istime))
        data(istime)=idft(data(istime));
    end
    if(any(isxy))
        data(isxy)=idft(data(isxy));
        data(isxy)=changeheader(data(isxy),'iftype','ixy');
    end
    
    % detail message
    if(verbose); disp('Finished Spectral Whitening of Record(s)'); end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
