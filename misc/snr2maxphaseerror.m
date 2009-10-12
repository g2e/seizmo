function [phase]=snr2maxphaseerror(snr)
%SNR2MAXPHASEERROR   Returns maximum narrow-band phase error based on SNR 
%
%    Usage:    phase=snr2maxphaseerror(snr)
%
%    Description: SNR2MAXPHASEERROR(SNR) calculates the maximum phase shift
%     that can be caused by noise when the signal to noise ratio is SNR.
%     The algorithm assumes the signal and noise are simple sinusoidal
%     functions with different phase and amplitude (basically limiting the
%     use of this function to narrowband applications).  It is also assumed
%     that the given SNR is the signal to noise ratio rather than signal
%     plus noise to noise ratio (on average S2NR & S+N2NR are equal, so
%     this is reasonable).  SNR must be a numeric array of positive reals.
%     Output is in radians!
%
%    Notes:
%     - Cycle skipping is another source of error not touched on here.  To
%       do this analysis, use the envelope and the snr to make an estimate
%       on the maximum number of cycles skipped.  It may be better to come
%       up with a method for handling cycle skipping rather than including
%       it in measurement error though.
%
%    Examples:
%     Get max effect of noise on P arrival phase:
%       Ptimes=getarrival(data,'P');
%       snr=quicksnr(data,Ptimes+[-100 -20],Ptimes+[-20 40])
%       period=10; % Assume 10sec is dominant period
%       arrivalerror=period/(2*pi).*snr2maxphaseerror(snr);
%
%    See also: QUICKSNR

%     Version History:
%        Aug. 22, 2009 - initial version
%        Sep.  8, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 06:55 GMT

% todo:

% make sure snr is positive real
if(~isreal(snr) || any(snr<0))
    error('seizmo:snr2maxphaseerror:badInput',...
        'SNR must be positive reals!');
end

% amplitude and phase
sz=size(snr);          % get array size for later
snr=snr(:);            % assure column vector
nsnr=numel(snr);       % save numel before expansion
theta=0:0.001*pi:pi;   % maybe a bit too dense...
snr=snr(:,ones(numel(theta),1)); % expand
theta=theta(ones(nsnr,1),:);     % expand
amp=sqrt(1+snr.^2+2.*snr.*cos(theta));
phase=atan2(sin(theta)./amp,(snr+cos(theta))./amp);
phase=reshape(max(abs(phase.')),sz); % get max and reshape back

end
