function [x]=rescale_snr(snr,inrange,outrange)
%RESCALE_SNR    Rescales SNR using an error function
%
%    Usage:    x=rescale_snr(snr)
%              x=rescale_snr(snr,inrange)
%              x=rescale_snr(snr,inrange,outrange)
%
%    Description: X=RESCALE_SNR(SNR) rescales SNR so it varies according to
%     an error function.  This is to de-emphasizing high signal-to-noise
%     data in inversions.  Output values in X range from 1 to 2.  Values
%     less than 1 in SNR will have a corresponding value near 1 in X while
%     values above 100 in SNR will correspond to values near 2 in X.  The
%     strongest transition in X occurs for SNR from 30 to 70 where X goes
%     from about 1.1 to 1.9.
%
%     X=RESCALE_SNR(SNR,INRANGE) sets the SNR range over which the X
%     output goes from 1 to 2.  The default is [1 100].  INRANGE must be a
%     1x2 vector of positive reals.
%
%     X=RESCALE_SNR(SNR,INRANGE,OUTRANGE) sets the output range of X to
%     OUTRANGE.  The default is [1 2].
%
%    Notes:
%     - This is mainly for using SNR values in weighting schemes.
%     - TTSOLVE makes use of RESCALE_SNR to avoid overemphasizing high-SNR
%       results.
%
%    Examples:
%     Plot up some random signal-to-noise ratio vs rescaled values:
%      snr=20*rand(100,1);
%      plot(snr,rescale_snr(snr),'x')
%
%    See also: TTSOLVE, SNR2MAXPHASEERROR

%     Version History:
%        Mar. 12, 2010 - initial version
%        Mar. 14, 2010 - doc update, limit to range
%        Sep. 13, 2010 - nargchk fix
%        Jan. 13, 2011 - readjusted input scaling to reflect signal to
%                        noise ratios between the rms of the noise and one
%                        half the peak 2 peak amplitude of the signal
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 13, 2011 at 01:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% default ranges
if(nargin<2 || isempty(inrange)); inrange=[1 100]; end
if(nargin<3 || isempty(outrange)); outrange=[1 2]; end

% check inputs
if(~isreal(snr) || any(snr(:)<0))
    error('seizmo:rescale_snr:badInput',...
        'SNR must be an array of positive reals!');
elseif(~isequal(size(inrange),[1 2]) || ~isreal(inrange) || any(inrange<0))
    error('seizmo:rescale_snr:badInput',...
        'INLIMIT must be a 1x2 positive real array!');
elseif(~isequal(size(outrange),[1 2]) || ~isnumeric(outrange))
    error('seizmo:rescale_snr:badInput',...
        'OUTLIMIT must be a 1x2 numeric array!');
end

% rescale/recenter snr
x=4/abs(diff(inrange))*(snr-mean(inrange));

% apply error function
x=erf(x);

% rescale/recenter x to outrange
x=x*abs(diff(outrange))/2+mean(outrange);

end
