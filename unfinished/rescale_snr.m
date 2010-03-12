function [x]=rescale_snr(snr,inlimit,outlimit)
%RESCALE_SNR    Rescales SNR using an error function
%
%    Usage:    x=rescale_snr(snr)
%              x=rescale_snr(snr,inlimit)
%              x=rescale_snr(snr,inlimit,outlimit)
%
%    Description: X=RESCALE_SNR(SNR) rescales SNR so it varies according to
%     an error function.  This is to de-emphasizing high signal-to-noise
%     data in inversions.  Output values in X range from 1 to 2.  Values
%     less than 3 in SNR will have a corresponding value near 1 in X while
%     values above 10 in SNR will correspond to values near 2 in X.  The
%     strongest transition in X occurs for SNR from 5 to 8 where X goes
%     from about 1.1 to 1.9.
%
%     X=RESCALE_SNR(SNR,INLIMIT) sets the SNR limits over which the X
%     output goes from 1 to 2.  The default is [3 10].  INLIMIT must be a
%     1x2 vector of positive reals.
%
%     X=RESCALE_SNR(SNR,INLIMIT,OUTLIMIT) sets the output range of X to
%     OUTLIMIT.  The default is [1 2].
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 12, 2010 at 16:05 GMT

% todo:

% check nargin
msg=nargchk(1,3,nargin);
if(~isempty(msg)); error(msg); end

% default limits
if(nargin<2 || isempty(inlimit)); inlimit=[3 10]; end
if(nargin<3 || isempty(outlimit)); outlimit=[1 2]; end

% check inputs
if(~isreal(snr) || any(snr(:)<0))
    error('seizmo:rescale_snr:badInput',...
        'SNR must be an array of positive reals!');
elseif(~isequal(size(inlimit),[1 2]) || ~isreal(inlimit) || any(inlimit<0))
    error('seizmo:rescale_snr:badInput',...
        'INLIMIT must be a 1x2 positive real array!');
elseif(~isequal(size(outlimit),[1 2]) || ~isnumeric(outlimit))
    error('seizmo:rescale_snr:badInput',...
        'OUTLIMIT must be a 1x2 numeric array!');
end

% rescale/recenter snr
x=4/abs(diff(inlimit))*(snr-mean(inlimit));

% apply error function
x=erf(x);

% rescale/recenter x to outlimit
x=x*abs(diff(outlimit))/2+mean(outlimit);

end
