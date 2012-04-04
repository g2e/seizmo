function [phase]=snr2phaseerror(snr,func)
%SNR2PHASEERROR   Returns narrow-band phase error based on SNR 
%
%    Usage:    phase=snr2phaseerror(snr)
%              phase=snr2phaseerror(snr,func)
%
%    Description:
%     SNR2PHASEERROR(SNR) calculates the phase shift that can be caused by
%     noise when the signal to noise ratio is SNR in a monochromatic case.
%     Please note that this is a 1 standard deviation measure of the phase
%     shifts (you can do other measures using the second usage form below).
%     The algorithm assumes the signal and noise are simple sinusoidal
%     functions of the same frequency with different phase and amplitude
%     (basically limiting the use of this function to narrowband
%     applications).  It is also assumed that the given SNR is the signal
%     to noise ratio rather than signal plus noise to noise ratio (on
%     average S/N & (S+N)/N are equal, so this is reasonable).  SNR must be
%     a numeric array of positive reals.  Output is in radians!
%
%     SNR2PHASEERROR(SNR,FUNC) alters the phase error function.  The
%     default is @std.  Other valid types are @max, @var, @median.  Note
%     that @median will probably always return 0.
%
%    Notes:
%     - Cycle skipping is another source of error not touched on here.  To
%       do this analysis, use the envelope and the snr to make an estimate
%       on the maximum number of cycles skipped.  It may be better to come
%       up with a method for handling cycle skipping rather than including
%       it in measurement error though.
%
%    Examples:
%     % Observed max travel-time error from noise on some 10s P waves:
%     P=findpicks(data,'P',1);
%     snr=quicksnr(data,P+[-100 -20],P+[-20 40])
%     period=10; % Assume 10sec is dominant period
%     maxtterror=period/(2*pi).*snr2phaseerror(snr,@max);
%
%    See also: QUICKSNR, TTSOLVE

%     Version History:
%        Aug. 22, 2009 - initial version
%        Sep.  8, 2009 - minor doc update
%        Jan. 18, 2011 - doc update
%        Feb. 12, 2011 - minor doc update, 2nd input, now std is used
%                        rather than max, name changed to snr2phaseerror
%        Mar. 15, 2012 - update example
%        Mar. 23, 2012 - simplified computation (removed amplitude term)
%                        with output changed only when snr=1 (now max=pi/2
%                        rather than max=pi) which is correct anyway
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 23, 2012 at 06:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% make sure snr is positive real
if(~isreal(snr) || any(snr<0))
    error('seizmo:snr2phaseerror:badInput',...
        'SNR must be positive reals!');
end

% default & check fun
if(nargin==1 || isempty(func)); func=@std; end
if(~isa(func,'function_handle'))
    error('seizmo:snr2phaseerror:badInput',...
        'FUNC must be a function handle!');
end

% amplitude and phase
sz=size(snr);          % get array size for later
snr=snr(:);            % assure column vector
nsnr=numel(snr);       % save numel before expansion
theta=0:0.1*pi:2*pi; % maybe a bit too dense...
snr=snr(:,ones(numel(theta),1)); % expand
theta=theta(ones(nsnr,1),:);     % expand
phase=atan2(sin(theta),snr+cos(theta));
phase=reshape(func(phase.'),sz); % get stddev and reshape back

end
