function [data]=whitensp(data,halfwindow,varargin)
%WHITENSP    Spectral whitening/normalization of SEIZMO spectral records
%
%    Usage:    data=whitensp(data)
%              data=whitensp(data,halfwindow)
%              data=whitensp(...,'optionname',optionvalue,...)
%
%    Description: WHITENSP(DATA) will perform spectral normalization (aka
%     whitening) on spectral records in the SEIZMO structure DATA.
%     Normalization is performed by dividing the complex spectrum by a
%     smoothed version of the amplitude spectrum.  Smoothing utilizes a 41
%     sample sliding mean.  This is NOT equivalent to the 'whiten' command
%     in SAC (see function PREWHITEN).  This operation is particularly
%     suited for ambient noise studies.  WHITESP differs from WHITEN in
%     that it expects and returns spectral records rather than time series.
%     Records are returned as Real-Imaginary.
%
%     WHITENSP(DATA,N) allows changing the half-width of the smoothing
%     window to 2N+1.  The default N is 20.  See function SLIDINGABSMEAN
%     for more information.
%
%     WHITENSP(...,'OPTIONNAME',OPTIONVALUE,...) will pass sliding average
%     options on to the SLIDINGABSMEAN call.  See SLIDINGABSMEAN for more
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
%      plotsp1([data(1) whitensp(data(1))],'am')
%
%    See also: whiten, slidingabsmean, prewhiten, unprewhiten

%     Version History:
%        June 11, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2009 at 16:00 GMT

% todo:

% check number of inputs
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% get some header fields
leven=getlgc(data,'leven');
iftype=getenumdesc(data,'iftype');

% require evenly-spaced time series, general x vs y
if(any(~strcmpi(leven,'true')))
    error('seizmo:whitensp:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Spectral File-Real/Imag')...
        & ~strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:whitensp:illegalOperation',...
        'Illegal operation on non-spectral record!')
end

% default centered window half size
if(nargin==1 || isempty(halfwindow)); halfwindow=2; end

% get amph and rlim type records
amph=rlim2amph(data);
data=amph2rlim(data);

% fake amph records as rlim (to get by dividerecords checks/fixes)
amph=changeheader(amph,'iftype','irlim');

% get smoothed amplitude records
amph=slidingabsmean(amph,halfwindow,varargin{:});

% copy amplitude over phase
amph=seizmofun(amph,@(x)x(:,[1:2:end; 1:2:end]));

% divide complex by smoothed amplitude
data=dividerecords(data,amph);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
