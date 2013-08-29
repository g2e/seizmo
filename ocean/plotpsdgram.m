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
%     for each frequency.  This is useful for inspecting subtle variations
%     in the spectra.  The default DSPEC is false (does not plot the
%     deviatoric spectrum).
%
%     PLOTPSDGRAM(AX,...) plots the spectrogram(s) in the axes given by AX.
%
%     AX=PLOTPSDGRAM(...) returns the handle of the axes plotted in.
%
%    Notes:
%
%    Examples:
%     % Read and plot the PSD for bouy 45001 from a
%     % series of storms in the middle of October 2011:
%     file=which('45001w2011101300-2011102300.txt');
%     psdgram=read_ndbc_swden(file);
%     plotpsdgram(psdgram);
%
%    See also: READ_NDBC_SWDEN, NOISE_PSDGRAM, CHKPSDGRAM

%     Version History:
%        Feb. 27, 2013 - initial version
%        Apr. 15, 2013 - added dtlim, flim, dspec inputs
%        Aug.  9, 2013 - simple updates to allow psdgram to be an array
%        Aug. 28, 2013 - fixed example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% strip off axes handle(s)
ax=[];
if(isnumeric(varargin{1}))
    ax=varargin{1};
    varargin=varargin(2:end);
    for i=1:numel(ax)
        if(~ishghandle(ax(i),'axes'))
            warning('seizmo:plotpsdgram:badAx',...
                'AX has an invalid handle!  Using default AX instead.');
            ax=[];
        end
    end
end

% check struct
psdgram=varargin{1};
error(chkpsdgram(psdgram));
npsd=numel(psdgram);

% require same number of psdgrams & axes
if(~isempty(ax) && numel(ax)~=npsd)
    warning('seizmo:plotpsdgram:badAx',...
        ['AX must have the same number of handles as\n' ...
        'elements in PSDGRAM!  Using default AX instead.']);
    ax=[];
end

% default & check optional inputs
if(nargin<2 || isempty(varargin{2}))
    step0=[];
else
    step0=varargin{2};
    if(~isscalar(step0) || ~isnumeric(step0) || ~isreal(step0) || step0<=0)
        error('seizmo:plotpsdgram:badInput',...
            'STEP must be a real-valued scalar in seconds!');
    elseif(fix(step0)~=step0)
        error('seizmo:plotpsdgram:badInput',...
            'STEP is not allowed to have sub-second accuracy!');
    end
end
if(nargin<3 || isempty(varargin{3}))
    dtlim0=[];
else
    dtlim0=varargin{3};
    if(~isequal(size(dtlim0),[1 2]) || ~isnumeric(dtlim0) ...
            || ~isreal(dtlim0))
        error('seizmo:plotpsdgram:badInput',...
            'DTLIM must be a 1x2 array [DTMIN DTMAX] in DATENUM format!');
    elseif(diff(dtlim0)<0)
        error('seizmo:plotpsdgram:badInput',...
            'DTLIM values are not in proper order!');
    end
end
if(nargin<4 || isempty(varargin{4}))
    flim0=[];
else
    flim0=varargin{4};
    if(~isequal(size(flim0),[1 2]) || ~isnumeric(flim0) || ~isreal(flim0))
        error('seizmo:plotpsdgram:badInput',...
            'FLIM must be a 1x2 real-valued array of [FMIN FMAX] in Hz!');
    elseif(any(flim0<0))
        error('seizmo:plotpsdgram:badInput',...
            'FLIM values must be positive!');
    elseif(diff(flim0)<0)
        error('seizmo:plotpsdgram:badInput',...
            'FLIM values not in proper order!');
    end
end
if(nargin<5 || isempty(varargin{5}))
    dspec=false;
else
    dspec=varargin{5};
    if(~isscalar(dspec) || ~islogical(dspec))
        error('seizmo:plotpsdgram:badInput',...
            'DSPEC must be TRUE or FALSE!');
    end
end

% default axes setup
if(isempty(ax))
    % new figure
    fh=figure;
    ncols=fix(sqrt(npsd));
    nrows=ceil(npsd/ncols);
    ax=makesubplots(nrows,ncols,1:npsd,'align','parent',fh);
end

% loop over each psdgram/axes
for i=1:npsd
    % defaults
    if(isempty(step0))
        step=round(min(diff(psdgram(i).time))*86400);
    else
        step=step0;
    end
    if(isempty(dtlim0))
        dtlim=[min(psdgram(i).time) max(psdgram(i).time)];
    else
        dtlim=dtlim0;
    end
    if(isempty(flim0))
        flim=[min(psdgram(i).freq) max(psdgram(i).freq)];
    else
        flim=flim0;
    end
    
    % apply limits
    dtok=psdgram(i).time>=dtlim(1) & psdgram(i).time<=dtlim(2);
    fok=psdgram(i).freq>=flim(1) & psdgram(i).freq<=flim(2);
    psdgram(i).time=psdgram(i).time(dtok);
    psdgram(i).freq=psdgram(i).freq(fok);
    psdgram(i).spectra=psdgram(i).spectra(dtok,fok);
    
    % add gaps, 0=>nan, convert to dB
    tall=0:step:(max(psdgram(i).time)-min(psdgram(i).time))*86400;
    [tf,loc]=ismember(...
        round((psdgram(i).time-min(psdgram(i).time))*86400),tall);
    if(any(~tf)); die; end
    spec=nan(numel(tall),numel(psdgram(i).freq));
    spec(loc,:)=psdgram(i).spectra;
    spec(spec==0)=nan;
    spec=10*log10(spec);
    if(dspec); spec=spec-repmat(nanmedian(spec),numel(tall),1); end
    
    % plot
    h=imagesc(tall/86400+min(psdgram(i).time),psdgram(i).freq,spec',...
        'parent',ax(i));
    set(ax(i),'color','none',...
        'xtick',floor(min(psdgram(i).time)):10:ceil(max(psdgram(i).time)));
    set(h,'alphadata',~isnan(spec'));
    datetick(ax(i),'x','mmmdd','keeplimits','keepticks');
    grid(ax(i),'on');
    ylabel(ax(i),'Freq (Hz)');
    cbh=colorbar('peer',ax(i));
    ylabel(cbh,['db (10*log_{10}(' psdgram(i).units '))']);
    title(ax(i),psdgram(i).name);
end

% output if desired
if(nargout); varargout={ax}; end

end
