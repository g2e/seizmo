function [varargout]=plotpsdgram(varargin)
%PLOTPSDGRAM    Plots power spectral density spectrogram struct
%
%    Usage:    plotpsdgram(psdgram)
%              plotpsdgram(psdgram,step)
%              plotpsdgram(psdgram,step,dtlim,flim)
%              plotpsdgram(psdgram,step,dtlim,flim,dspec)
%              plotpsdgram(ax,...)
%              ax=plotpsdgram(...)
%
%    Description:
%     PLOTPSDGRAM(PSDGRAM) makes a spectrogram plot using the power
%     spectral density (PSD) info in struct PSDGRAM.  PSDGRAM should be as
%     output from READ_NDBC_SWDEN.  The time step is determined by the
%     minimum time interval between PSD points.
%
%     PLOTPSDGRAM(PSDGRAM,STEP) uses STEP as the time interval for the
%     spectrogram.  STEP is in seconds.  Sub-second accuracy is not
%     currently allowed.  A STEP less than the interval between PSD time
%     points will insert columns of NaNs which will appear gray in the
%     plot.
%
%     PLOTPSDGRAM(PSDGRAM,STEP,DTLIM,FLIM) limits the plotted spectrum to
%     the date-time range in DTLIM and the frequency range in FLIM.  DTLIM
%     must be in datenum format (see DATENUM, DATESTR, DATEVEC, etc) as
%     [DTMIN DTMAX].  FLIM must be in Hz as [FMIN FMAX].
%
%     PLOTPSDGRAM(PSDGRAM,STEP,DTLIM,FLIM,DSPEC) indicates whether or not
%     the plotted spectrum should be the deviation from the median value
%     for each frequency.  This is useful for looking for inspecting subtle
%     variations in the spectra.  The default DSPEC is false (do not plot
%     the deviatoric spectrum).
%
%     PLOTPSDGRAM(AX,...) plots the spectrogram in the axes AX.
%
%     AX=PLOTPSDGRAM(...) returns the handle of the axes plotted in.
%
%    Notes:
%
%    Examples:
%     % Bouy 45001 (Middle Lake Superior) spectrogram for 2011:
%     path=fileparts(which('read_ndbc_swden'));
%     s=read_ndbc_swden([path filesep '45001w2011.txt']);
%     plotpsdgram(s);
%
%    See also: READ_NDBC_SWDEN, NOISE_PSDGRAM, CHKPSDGRAM

%     Version History:
%        Feb. 27, 2013 - initial version
%        Apr. 15, 2013 - added dtlim, flim, dspec inputs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 15, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% strip off axes handle
[ax,varargin]=axparse(varargin{:});
if(isempty(ax)); ax=gca; end

% check struct
psdgram=varargin{1};
error(chkpsdgram(psdgram));

% defaults
if(nargin<2 || isempty(varargin{2}))
    step=round(min(diff(psdgram.time))*86400);
else
    step=varargin{2};
end
if(nargin<3 || isempty(varargin{3}))
    dtlim=[min(psdgram.time) max(psdgram.time)];
else
    dtlim=varargin{3};
end
if(nargin<4 || isempty(varargin{4}))
    flim=[min(psdgram.freq) max(psdgram.freq)];
else
    flim=varargin{4};
end
if(nargin<5 || isempty(varargin{5}))
    dspec=false;
else
    dspec=varargin{5};
end

% check optional inputs
if(~isscalar(step) || ~isnumeric(step) || ~isreal(step) || step<=0)
    error('seizmo:plotpsdgram:badInput',...
        'STEP must be a real-valued scalar in seconds!');
elseif(fix(step)~=step)
    error('seizmo:plotpsdgram:badInput',...
        'STEP is not allowed to have sub-second accuracy!');
elseif(~isequal(size(dtlim),[1 2]) || ~isnumeric(dtlim) || ~isreal(dtlim))
    error('seizmo:plotpsdgram:badInput',...
        'DTLIM must be a 1x2 array of [DTMIN DTMAX] in DATENUM format!');
elseif(diff(dtlim)<0)
    error('seizmo:plotpsdgram:badInput',...
        'DTLIM values are not in proper order!');
elseif(~isequal(size(flim),[1 2]) || ~isnumeric(flim) || ~isreal(flim))
    error('seizmo:plotpsdgram:badInput',...
        'FLIM must be a 1x2 real-valued array of [FMIN FMAX] in Hz!');
elseif(any(flim<0))
    error('seizmo:plotpsdgram:badInput',...
        'FLIM values must be positive!');
elseif(diff(flim)<0)
    error('seizmo:plotpsdgram:badInput',...
        'FLIM values not in proper order!');
elseif(~isscalar(dspec) || ~islogical(dspec))
    error('seizmo:plotpsdgram:badInput',...
        'DSPEC must be TRUE or FALSE!');
end

% apply limits
dtok=psdgram.time>=dtlim(1) & psdgram.time<=dtlim(2);
fok=psdgram.freq>=flim(1) & psdgram.freq<=flim(2);
psdgram.time=psdgram.time(dtok);
psdgram.freq=psdgram.freq(fok);
psdgram.spectra=psdgram.spectra(dtok,fok);

% add gaps, 0=>nan, convert to dB
tall=0:step:(max(psdgram.time)-min(psdgram.time))*86400;
[tf,loc]=ismember(round((psdgram.time-min(psdgram.time))*86400),tall);
if(any(~tf)); die; end
spec=nan(numel(tall),numel(psdgram.freq));
spec(loc,:)=psdgram.spectra;
spec(spec==0)=nan;
spec=10*log10(spec);
if(dspec); spec=spec-repmat(nanmedian(spec),numel(tall),1); end

% plot
h=imagesc(tall/86400+min(psdgram.time),psdgram.freq,spec','parent',ax);
set(ax,'color','none',...
    'xtick',floor(min(psdgram.time)):10:ceil(max(psdgram.time)));
set(h,'alphadata',~isnan(spec'));
datetick(ax,'x','mmmdd','keeplimits','keepticks');
grid(ax,'on');
ylabel(ax,'Freq (Hz)');
cbh=colorbar('peer',ax);
ylabel(cbh,['db (10*log_{10}(' psdgram.units '))']);
title(psdgram.name);

% output if desired
if(nargout); varargout={ax}; end

end
