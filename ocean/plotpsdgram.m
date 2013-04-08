function [varargout]=plotpsdgram(varargin)
%PLOTPSDGRAM    Plots power spectral density spectrogram struct
%
%    Usage:    plotpsdgram(psdgram)
%              plotpsdgram(psdgram,step)
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
%     currently allowed.
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
%    See also: READ_NDBC_SWDEN

%     Version History:
%        Feb. 27, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% strip off axes handle
[ax,varargin]=axparse(varargin{:});
if(isempty(ax)); ax=gca; end

% check struct
psdgram=varargin{1};
error(chkpsdgram(psdgram));

% default step
if(nargin<2 || isempty(varargin{2}))
    step=round(min(diff(psdgram.time))*86400);
else
    step=varargin{2};
end

% check step
if(~isscalar(step) || ~isnumeric(step) || ~isreal(step) || step<=0)
    error('seizmo:plotpsdgram:badStep',...
        'STEP must be a real-valued scalar in seconds!');
elseif(fix(step)~=step)
    error('seizmo:plotpsdgram:badStep',...
        'STEP is not allowed to have sub-second accuracy!');
end

% add gaps, 0=>nan, convert to dB
tall=0:step:(max(psdgram.time)-min(psdgram.time))*86400;
[tf,loc]=ismember(round((psdgram.time-min(psdgram.time))*86400),tall);
if(any(~tf)); die; end
spec=nan(numel(tall),numel(psdgram.freq));
spec(loc,:)=psdgram.spectra;
spec(spec==0)=nan;
spec=10*log10(spec);
%spec=spec-repmat(nanmedian(spec),numel(tall),1);

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

% output if desired
if(nargout); varargout={ax}; end

end
