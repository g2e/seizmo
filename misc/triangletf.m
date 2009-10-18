function [x]=triangletf(t,t0,hwidth,amp)
%TRIANGLETF    Returns a triangle time function
%
%    Usage:    x=triangletf(t,t0,hwidth)
%              x=triangletf(t,t0,hwidth,amp)
%
%    Description: X=TRIANGLETF(T,T0,HWIDTH) creates a triangle time
%     function centered at T0 (with an amplitude of 1) that is 0 at times
%     T0-HWIDTH and T0+HWIDTH.  Values outside of this range are 0.
%
%     X=TRIANGLETF(T,T0,HWIDTH,AMP) sets the peak amplitude of the triangle
%     (at time T0) to AMP.  The default value of AMP is 1.
%
%    Notes:
%
%    Examples:
%     Compare gaussian and triangle time functions
%     for convolving with synthetic seismic data:
%      plot(-15:0.1:15,gaussiantf(-15:0.1:15,0,10,...
%          exp(0.5)*sqrt(2),exp(0.5)/10/sqrt(pi)))
%      hold on
%      plot(-15:0.1:15,triangletf(-15:0.1:15,0,10,1/10))
%
%    See also: GAUSSIANTF, TAPERFUN, TRIANG

%     Version History:
%        Oct. 17, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2009 at 21:15 GMT

% todo:

% check nargin
msg=nargchk(2,4,nargin);
if(~isempty(msg)); error(msg); end

% defaults
if(nargin<4 || isempty(amp)); amp=1; end

% check inputs
if(~isnumeric(t))
    error('seizmo:gaussiantf:badInput','T must be a numeric array!');
elseif(~isscalar(t0) || ~isnumeric(t0))
    error('seizmo:gaussiantf:badInput','TO must be a numeric scalar!');
elseif(~isscalar(hwidth) || ~isnumeric(hwidth))
    error('seizmo:gaussiantf:badInput','HWIDTH must be a numeric scalar!');
elseif(~isscalar(amp) || ~isnumeric(amp))
    error('seizmo:gaussiantf:badInput','AMP must be a numeric scalar!');
end

% get triangle values
x=amp.*(1-abs((t-t0)./hwidth));
x(x<0)=0;

end
