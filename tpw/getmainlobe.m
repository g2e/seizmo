function [f,a]=getmainlobe(f0,fs,swin,tprfrac,zpad,show)
%GETMAINLOBE    Returns the main lobe of a windowed sinusoid
%
%    Usage:    [f,a]=getmainlobe(f0,fs,swin,tprfrac,zpad)
%              [f,a]=getmainlobe(f0,fs,swin,tprfrac,zpad,show)
%
%    Description:
%     [F,A]=GETMAINLOBE(F0,FS,SWIN,TPRFRAC,ZPAD) extracts the frequency and
%     amplitude of the main lobe of a windowed sinusoid.  The sinusoid's
%     frequency is given by F0 (in Hz) and has unit amplitude.  The
%     sampling frequency is given by FS (in Hz).  FS must be greater than
%     twice the sinusoid's frequency to avoid aliasing.  SWIN gives the
%     boxcar window limits ([START END] in seconds), TPRFRAC gives the
%     cosine taper halfwidth as a fraction of the boxcar window size, and
%     ZPAD gives amouunt of zero-padding before and after the window
%     ([BEFORE AFTER] in seconds).  Given these specifications, the
%     sinusoid is boxcar windowed, cosine tapered on both ends, and zero
%     padded.  Finally, the windowed sinusoid is converted to the frequency
%     domain and the lobe with the strongest peak amplitude is extracted.
%
%     [F,A]=GETMAINLOBE(F0,FS,SWIN,TPRFRAC,ZPAD,SHOW) toggles the plotting
%     of the sinusoid using SHOW.  Plots are of the sinusoid after
%     windowing, tapering and padding and of the corresponding spectral
%     amplitudes.  SHOW must be TRUE (make plots) or FALSE (no plotting,
%     the default).
%
%    Notes:
%     - sinusoid is a sine function w/ no phase offset (amp. at 0s is 0)
%
%    Examples:
%     % Get the main lobe for a 100s period sinusoid sampled at 1Hz with a
%     % 1000s window with 200 seconds of cosine tapering on each end and
%     % 1000s of zero padding on both ends:
%     [f,a]=getmainlobe(1/100,1,[1000 2000],200/1000,[1000 1000]);
%
%     % Compare that with a non-tapered version:
%     [f1,a1]=getmainlobe(1/100,1,[1000 2000],0,[1000 1000]);
%     plot(f,a,f1,a1)
%
%    See also: RAYLEIGH_2D_PLANE_WAVE_KERNELS, SMOOTH2D, READKERNELS,
%              WRITEKERNELS, MAKEKERNELS, PLOTKERNELS

%     Version History:
%        Feb.  4, 2010 - initial version
%        July  9, 2010 - fixed see also section, fixed nargchk
%        Mar. 24, 2012 - minor doc update, plot calls use handles now
%        May  30, 2012 - pow2pad=0 by default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2012 at 17:10 GMT

% todo:

% check nargin
error(nargchk(5,6,nargin));

% check inputs
if(nargin==5 || isempty(show)); show=false; end
if(~isscalar(f0) || ~isreal(f0) || f0<=0)
    error('seizmo:getmainlobe:badInput',...
        'F0 must be a positive real scalar (in Hz)!');
elseif(~isscalar(fs) || ~isreal(fs) || fs<=2*f0)
    error('seizmo:getmainlobe:badInput',...
        'FS must be a positive real scalar (in Hz) & FS>(2*F0) !');
elseif(numel(swin)~=2 || ~isreal(swin))
    error('seizmo:getmainlobe:badInput',...
        'SWIN must be an 2-element array of [START END] (in sec)!');
elseif(~isscalar(tprfrac) || ~isreal(tprfrac) || tprfrac<0 || tprfrac>1)
    error('seizmo:getmainlobe:badInput',...
        'TPRFRAC must be a positive real scalar from 0 to 1 !');
elseif(numel(zpad)~=2 || ~isreal(zpad) || any(zpad<0))
    error('seizmo:getmainlobe:badInput',...
        'ZPAD must be an 2-element array of [PREPAD POSTPAD] (in sec)!');
elseif(~isscalar(show) || ~islogical(show))
    error('seizmo:getmainlobe:badInput',...
        'SHOW must be a scalar logical!');
end

% sort window
swin=sort(swin);

% generate boxcar windowed sine function
t=swin(1):1/fs:swin(2);
npts=numel(t);
x=sin(2*pi*f0*t);
%figure; plot(t,x);

% cosine taper the window
ntpts=round(npts*tprfrac);
tpr=taperfun('cosine',(1:ntpts)/ntpts,[0 1]);
x(1:ntpts)=x(1:ntpts).*tpr;
x(end-ntpts+1:end)=x(end-ntpts+1:end).*tpr(end:-1:1);
%figure; plot(t,x);

% zero pad
nzp=round(zpad*fs);
t=[t(1)-(nzp(1):-1:1)/fs t t(end)+(1:nzp(2))/fs];
x=[zeros(1,nzp(1)) x zeros(1,nzp(2))];
npts=nzp(1)+nzp(2)+npts;
if(show)
    fh(1)=figure;
    ax(1)=axes('parent',fh(1));
    plot(ax(1),t,x);
    xlabel(ax(1),'TIME (s)');
    ylabel(ax(1),'NON-DIMENSIONAL AMPLITUDE');
    title(ax(1),['WINDOWED SINUSOID - ' sprintf('PERIOD: %gS',1/f0)]);
    grid(ax(1),'on');
end

% fft
nfft=2^nextpow2(npts);
c=fft(x,nfft)/fs;

% get frequency & amplitudes
f=fs/nfft*(0:nfft/2);
a=2*abs(c(1:nfft/2+1));
if(show)
    fh(2)=figure;
    ax(2)=axes('parent',fh(2));
    loglog(ax(2),f,a);
    xlabel(ax(2),'FREQUENCY (Hz)');
    ylabel(ax(2),'NON-DIMENSIONAL SPECTRAL AMPLITUDE');
    title(ax(2),['WINDOWED SINUSOID SPECTRA - ' ...
        sprintf('PERIOD: %gS',1/f0)]);
    grid(ax(2),'on');
end

% extract the main lobe
[i,i]=max(a);
da=diff(a);
i1=find(da(1:i-1)<0,1,'last')+1;
i2=i+find(da(i+1:end)>0,1,'first');
if(isempty(i1)); i1=1; end
if(isempty(i2)); i2=numel(a); end
f=f(i1:i2);
a=a(i1:i2);
%figure; plot(f,a);

end
