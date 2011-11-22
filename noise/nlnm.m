function [db,freq]=nlnm(freq,units)
%NLNM    Returns the power spectral density of the USGS low noise model
%
%    Usage:    [db,freq]=nlnm
%              db=nlnm(freq)
%              db=nlnm(freq,units)
%
%    Description:
%     [DB,FREQ]=NLNM returns the Peterson 1993 New Low Noise Model (NLNM)
%     for the Earth's seismic background from frequencies of 1e-5 to 10 Hz.
%     The noise model is in units of decibels relative to (1 m/s^2)^2/Hz.
%
%     DB=NLNM(FREQ) returns the NLNM at the frequencies specified.  Note
%     they must be in the range 1e-5 to 10 Hz.
%
%     DB=NLNM(FREQ,UNITS) allows getting the NLNM in acceleration,
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
%    See also: NHNM

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
    error('seizmo:nlnm:badInput',...
        'FREQ must be real-valued and be in Hz between .00001 & 10 Hz!');
elseif(~ischar(units) || size(units,1)~=1 ...
        || ~any(strcmpi(units,valid.units)))
    error('seizmo:nlnm:badInput',...
        'UNITS must be ''acc'', ''vel'', or ''dis''!');
end

% NLNM
pab=[.1     -162.36 5.64
     .17    -166.7  0
     .4     -170    -8.3
     .8     -166.4  28.9
     1.24   -168.6  52.48
     2.4    -159.98 29.81
     4.3    -141.1  0
     5      -71.36  -99.77
     6      -97.26  -66.49
     10     -132.18 -31.57
     12     -205.27 36.16
     15.6   -37.65  -104.33
     21.9   -114.37 -47.1
     31.6   -160.58 -16.28
     45     -187.5  0
     70     -216.47 15.7
     101    -185    0
     154    -168.34 -7.61
     328    -217.43 11.9
     600    -258.28 26.6
     10000  -346.88 48.75];
 
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
