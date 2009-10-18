function [x]=gaussiantf(t,t0,hwidth,stddev,amp)
%GAUSSIANTF    Returns a gaussian time function
%
%    Usage:    x=gaussiantf(t,t0,hwidth)
%              x=gaussiantf(t,t0,hwidth,stddev)
%              x=gaussiantf(t,t0,hwidth,stddev,amp)
%
%    Description: X=GAUSSIANTF(T,T0,HWIDTH) samples a gaussian curve
%     centered at time TO at the times in array T.  The gaussian curve is
%     further defined to falloff to a value of 1/e at a time of HWIDTH from
%     TO.  Both TO and HWIDTH must be numeric scalars.  T must be a numeric
%     array and X is an equal sized array containing the associated
%     gaussian values.
%
%     X=GAUSSIANTF(T,T0,HWIDTH,STDDEV) sets the falloff value STDDEV.  The
%     default value is sqrt(2), which causes the gaussian to drop off to
%     1/e at a time distance of HWIDTH from T0.  See the Notes section for
%     the gaussian formula.
%
%     X=GAUSSIANTF(T,T0,HWIDTH,STDDEV,AMP) sets the peak amplitude of the
%     gaussian curve at T0.  The default value of AMP is 1.
%
%    Notes:
%     - Formula for returned gaussian values:
%
%                   1     /           T-T0  \  2
%                - ___ * | STDDEV * ________ |
%                   2     \          HWIDTH /
%     X = AMP * e
%
%    Examples:
%     Compare gaussian and triangle time functions
%     for convolving with synthetic seismic data:
%      plot(-15:0.1:15,gaussiantf(-15:0.1:15,0,10,...
%          exp(0.5)*sqrt(2),exp(0.5)/10/sqrt(pi)))
%      hold on
%      plot(-15:0.1:15,triangletf(-15:0.1:15,0,10,1/10))
%
%    See also: TRIANGLETF, TAPERFUN, GAUSSWIN

%     Version History:
%        Oct. 11, 2009 - initial version
%        Oct. 17, 2009 - new example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2009 at 21:30 GMT

% todo:

% check nargin
msg=nargchk(3,5,nargin);
if(~isempty(msg)); error(msg); end

% defaults
if(nargin<4 || isempty(stddev)); stddev=sqrt(2); end
if(nargin<5 || isempty(amp)); amp=1; end

% check inputs
if(~isnumeric(t))
    error('seizmo:gaussiantf:badInput','T must be a numeric array!');
elseif(~isscalar(t0) || ~isnumeric(t0))
    error('seizmo:gaussiantf:badInput','TO must be a numeric scalar!');
elseif(~isscalar(hwidth) || ~isnumeric(hwidth))
    error('seizmo:gaussiantf:badInput','HWIDTH must be a numeric scalar!');
elseif(~isscalar(stddev) || ~isnumeric(stddev))
    error('seizmo:gaussiantf:badInput','STDDEV must be a numeric scalar!');
elseif(~isscalar(amp) || ~isnumeric(amp))
    error('seizmo:gaussiantf:badInput','AMP must be a numeric scalar!');
end

% get gaussian values
x=amp.*exp(-1./2.*(stddev.*(t-t0)./hwidth).^2);

end
