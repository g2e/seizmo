function [x]=scaled_erf(x,inrange,outrange)
%SCALED_ERF    Scale input using an error function
%
%    Usage:    y=scaled_erf(x)
%              y=scaled_erf(x,inrange)
%              y=scaled_erf(x,inrange,outrange)
%
%    Description:
%     Y=SCALED_ERF(X) rescales X so it varies according to an error
%     function.  Output values in X range from -1 to 1.  X values less than
%     -2 will have a corresponding Y value near -1 while X values above 2
%     will correspond to Y values near 1.  The strongest transition in Y
%     occurs for X values from -1 to 1 where Y goes from about -0.8 to 0.8.
%     This is no different from Y=ERF(X).
%
%     Y=SCALED_ERF(X,INRANGE) sets the X range over which the Y output goes
%     from -1 to 1.  The default is [-2 2].  INRANGE must be a 1x2 vector
%     of reals.
%
%     Y=SCALED_ERF(X,INRANGE,OUTRANGE) sets the output range of Y to
%     OUTRANGE.  The default is [-1 1].
%
%    Notes:
%
%    Examples:
%     % Scale some random signal-to-noise ratios
%     % (varying between 1 & 20) to between 1 & 2:
%     snr=19*rand(100,1)+1;
%     plot(snr,scaled_erf(snr,[1 20],[1 2]),'x')
%
%    See also: ERF, SNR2MAXPHASEERROR

%     Version History:
%        Mar. 12, 2010 - initial version
%        Mar. 14, 2010 - doc update, limit to range
%        Sep. 13, 2010 - nargchk fix
%        Jan. 13, 2011 - readjusted input scaling to reflect signal to
%                        noise ratios between the rms of the noise and one
%                        half the peak 2 peak amplitude of the signal
%        Jan. 18, 2011 - no longer used by ttsolve, renamed from
%                        rescaled_snr to scaled_erf, changed ranges
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 01:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% default ranges
if(nargin<2 || isempty(inrange)); inrange=[-2 2]; end
if(nargin<3 || isempty(outrange)); outrange=[erf(-2) erf(2)]; end

% check inputs
if(~isreal(x) || any(x(:)<0))
    error('seizmo:scaled_erf:badInput',...
        'X must be an array of positive reals!');
elseif(~isequal(size(inrange),[1 2]) || ~isreal(inrange))
    error('seizmo:scaled_erf:badInput',...
        'INLIMIT must be a 1x2 real array!');
elseif(~isequal(size(outrange),[1 2]) || ~isnumeric(outrange))
    error('seizmo:scaled_erf:badInput',...
        'OUTLIMIT must be a 1x2 numeric array!');
end

% rescale/recenter x
x=4/abs(diff(inrange))*(x-mean(inrange));

% apply error function
x=erf(x);

% rescale/recenter x to outrange
x=x*abs(diff(outrange))/2+mean(outrange);

end
