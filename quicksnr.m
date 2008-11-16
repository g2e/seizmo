function [snr]=qcksnr(data,reftimes,noswin,sigwin)
%QCKSNR    Quick estimation of SNR for SEIZMO data records
%
%    Description: Simply estimates the signal to noise ratio for SEIZMO 
%     data records by calculating the ratio of the maximum amplitudes of
%     two data windows that represent the 'noise' and the 'signal'.  Timing
%     parameters all need to be numeric - more advanced windowing should be
%     done explicitly with cutim.
%
%      reftimes - vector of reference times (eg phase arrival times)
%      noswin - [ws we] - noise window relative start and end times
%      sigwin - [ws we] - signal window relative start and end times
%
%    Usage: [snr]=qcksnr(data,reftimes,noisewindow,signalwindow)
%
%    Examples:
%     To get SNR estimates of P (assuming times are stored in header):
%       Ptimes=pullarr(data,'P')
%       snr=qcksnr(data,Ptimes,[-100 -20],[-20 40])
%
%    See also: pullarr, cutim, gnrm, nrm

% check nargin
error(nargchk(4,4,nargin))

% check data structure
error(seischk(data,'x'))

% get ratio of max amplitudes of windows
snr=gnrm(cutim(data,'z',reftimes+sigwin(:,1),'z',reftimes+sigwin(:,2))) ./ ...
gnrm(cutim(data,'z',reftimes+noswin(:,1),'z',reftimes+noswin(:,2)));
snr=snr(:); % make column vector

end
