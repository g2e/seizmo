function [db,freq]=nhnm(freq,units)
%NHNM    Returns the power spectral density of the USGS high noise model
%
%    Usage:    [db,freq]=nhnm
%              db=nhnm(freq)
%              db=nhnm(freq,units)
%
%    Description:
%     [DB,FREQ]=NHNM returns the Peterson 1993 New High Noise Model (NHNM)
%     for the Earth's seismic background from frequencies of 1e-5 to 10 Hz.
%     The noise model is in units of decibels relative to (1 m/s^2)^2/Hz.
%
%     DB=NHNM(FREQ) returns the NHNM at the frequencies specified.  Note
%     they must be in the range 1e-5 to 10 Hz.
%
%     DB=NHNM(FREQ,UNITS) allows getting the NHNM in acceleration,
%     velocity, or displacement.  The default UNITS is acceleration.  UNITS
%     must be a string such as 'acc', 'vel' or 'dis'.
%
%    Notes:
%
%    Examples:
%     % Make a plot showing the NLNM & NHNM
%     % over their entire frequency range:
%     fh=figure;
%     ax=axes('parent',fh);
%     plot(ax,10.^(-5:.1:1),nlnm,'g.-',10.^(-5:.1:1),nhnm,'r.-',...
%         'linewidth',2);
%     grid(ax,'on');
%     set(ax,'xscale','log');
%     xlabel(ax,'Frequency (Hz)');
%     ylabel(ax,'Power Spectral Density, dB (rel. 1 m^2/s^4/Hz)');
%     title(ax,'Peterson 1993 Global Background Noise Models');
%     legend(ax,'NLNM','NHNM');
%
%    See also: NLNM

%     Version History:
%        Nov.  9, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  9, 2011 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% defaults (100000s to 10Hz in acceleration)
if(nargin<1 || isempty(freq)); freq=10.^(-5:.1:1); end
if(nargin<2 || isempty(units)); units='acc'; end

% check inputs
valid.units={'acc' 'vel' 'dis' 'a' 'v' 'd' 'accel' 'velo' 'disp' ...
    'acceleration' 'velocity' 'displacement' 'none'};
if(~isnumeric(freq) || ~isreal(freq) ...
        || any(freq(:)<((1e-5)-eps) | freq(:)>10))
    error('seizmo:nhnm:badInput',...
        'FREQ must be real-valued and be in Hz between .00001 & 10 Hz!');
elseif(~ischar(units) || size(units,1)~=1 ...
        || ~any(strcmpi(units,valid.units)))
    error('seizmo:nhnm:badInput',...
        'UNITS must be ''acc'', ''vel'', or ''dis''!');
end

% NHNM
pab=[.1     -108.73 -17.23
     .22    -150.34 -80.5
     .32    -122.31 -23.87
     .8     -116.85 32.51
     3.8    -108.48 18.08
     4.6    -74.66  -32.95
     6.3    .66     -127.18
     7.9    -93.37  -22.42
     15.4   73.54   -162.98
     20     -151.52 10.01
     354.8  -206.66 31.63];
 
% frequency to period
period=1./freq;

% get indices of frequencies
% - just do it loop style cause that is simpler
sz=size(period);
idx=nan(sz);
for i=1:numel(period)
    idx(i)=sum(period(i)>=pab(:,1));
end

% get psd
db=reshape(pab(idx,2),sz)+reshape(pab(idx,3),sz).*log10(period);
switch units
    case {'v' 'vel' 'velo' 'velocity'}
        db=db+20*log10(period/(2*pi));
    case {'d' 'dis' 'disp' 'displacement' 'none'}
        db=db+20*log10((period/(2*pi)).^2);
end

end
