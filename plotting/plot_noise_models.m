function [varargout]=plot_noise_models(varargin)
%PLOT_NOISE_MODELS    Plots noise models (NLNM & NHNM)
%
%    Usage:    plot_noise_models
%              plot_noise_models(units)
%              plot_noise_models(ax,...)
%              h=plot_noise_models(...)
%
%    Description:
%     PLOT_NOISE_MODELS plots the Peterson 1993 New Low Noise Model (NLNM)
%     and New High Noise Model (NHNM) for the Earth's seismic background
%     from frequencies of 1e-5 to 10 Hz.  The noise model is in units of
%     decibels relative to (1 m/s^2)^2/Hz (which corresponds to the power
%     spectral density for acceleration seismograms).  Note you will have
%     to convert seismograms from nm/s^2 to m/s^2 (divide by 1e9)!
%
%     PLOT_NOISE_MODELS(UNITS)  allows getting the NLNM in acceleration,
%     velocity, or displacement.  The default UNITS is acceleration.  UNITS
%     must be a string such as 'acc', 'vel' or 'dis'.
%
%     PLOT_NOISE_MODELS(AX,...) indicates which plot(s) to draw the noise
%     models in.
%
%     H=PLOT_NOISE_MODELS(...) returns the handles of the noise model lines
%     as a Nx2 array where N is the number of axes drawn in.  H(:,1) gives
%     the low noise model handles and H(:,2) gives the high noise model
%     handles.
%
%    Notes:
%     - Does not work with PLOTSPECTRA0 or SPECTRASECTION
%     - Only valid if drawn against power spectra (use 'CMP','PW' option in
%       PLOTSPECTRA1/2) and the units are correct.  See the Examples
%       section below
%
%    Examples:
%     % Read in some data, remove response to acceleration, convert to
%     % m/s^2 from nm/s^2, convert to the frequency domain and plot power
%     % spectra against noise models:
%     data=readseizmo('some/files*');
%     data=removesacpz(data,'units','acc');
%     data=divide(data,1e9); % nm/s^2 to m/s^2
%     data=dft(data);
%     ax=plotspectra2(data,'cmp','pw');
%     plot_noise_models(ax);
%
%    See also: NLNM, NHNM, POWERSPECTRALDENSITY, KEEPPW, PLOTSPECTRA1,
%              PLOTSPECTRA2

%     Version History:
%        Feb.  6, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  6, 2012 at 23:00 GMT

% todo:

% get axes
if(nargin && isnumeric(varargin{1}))
    ax=ishghandle(varargin{1},'axes');
    if(any(ax))
        ax=varargin{1}(ax);
    else % not handles...
        error('seizmo:plot_noise_models:badInput',...
            'AX did not contain any valid axes handles!');
    end
    varargin=varargin(2:end);
    nargs=nargin-1;
else
    ax=gca; % use current
    nargs=nargin;
end

% check nargin
error(nargchk(0,1,nargs));

% check units
if(~nargs || isempty(varargin{1}))
    units='acc';
else
    units=varargin{1};
end

% draw on each axis
[dbl,fl]=nlnm([],units);
[dbh,fh]=nhnm([],units);
h=nan(numel(ax),2);
for i=1:numel(ax)
    held=ishold(ax(i));
    if(~held); hold(ax(i),'on'); end
    h(i,1)=plot(ax(i),fl,dbl,'color',[.5 .7 .5],'linewidth',2,...
        'displayname','NLNM','tag','NLNM');
    h(i,2)=plot(ax(i),fh,dbh,'color',[.7 .5 .5],'linewidth',2,...
        'displayname','NHNM','tag','NHNM');
    if(~held); hold(ax(i),'off'); end
end

if(nargout); varargout{1}=h; end

end
