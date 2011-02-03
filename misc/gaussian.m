function [p]=gaussian(x,mu,sigma)
%GAUSSIAN    Returns values from a Gaussian PDF
%
%    Usage:    p=gaussian(x)
%              p=gaussian(x,mu,sigma)
%
%    Description:
%     P=GAUSSIAN(X) returns the values of the standard normal (Gaussian)
%     distribution with zero mean and a stardard deviation of 1.  The
%     standard normal curve is a continuous probability density function
%     useful for quantifying and analyzing measurements with random error.
%     The probability density function is:
%                                          2
%                                 / -(X-MU)  \
%                                | __________ |
%           /                 \  |         2  |
%          |         1         |  \ 2*SIGMA  /
%      P = | _________________ | e
%          |          ______   |
%           \  SIGMA*V 2*PI   /
%
%     where MU is the mean of the random distribution (identifiable as the
%     peak location in the distribution and SIGMA is the standard deviation
%     and is a measure of the width of the distribution (identifiable as
%     the distance of the inflection points in the curve from the mean).
%     For the normal standard distribution MU=0 & SIGMA=1.
%
%     P=GAUSSIAN(X,MU,SIGMA) allows adjusting the mean MU and standard
%     deviation SIGMA of the gaussian distribution.  X, MU, & SIGMA must be
%     equal sized or scalars.
%
%    Notes:
%     - The area under the gaussian curve is equal to 1.  This is quite
%       useful in convolution as it preserves the energy in the process.
%
%    Examples:
%     % Create and plot several gaussian distributions:
%     x=(-5:.02:5)';
%     p0=gaussian(x);
%     p1=gaussian(x,-1);
%     p2=gaussian(x,0,2);
%     p3=gaussian(x,0,.5);
%     fh=figure; ax=axes('parent',fh);
%     plot(ax,x,[p0 p1 p2 p3],'linewidth',2);
%     xlabel(ax,'X');
%     ylabel(ax,'Probability Density');
%     legend(ax,{'\mu=0  \sigma=1' '\mu=-1 \sigma=1' ...
%                '\mu=0  \sigma=2' '\mu=0  \sigma=0.5'})
%     grid(ax,'on');
%
%    See also: GAUSSIANTF, TRIANGLETF

%     Version History:
%        Feb.  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  1, 2011 at 21:30 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% defaults
if(nargin<2 || isempty(mu)); mu=0; end
if(nargin<3 || isempty(sigma)); sigma=1; end

% check inputs
if(~isnumeric(x) || ~isnumeric(mu) || ~isnumeric(sigma))
    error('seizmo:gaussian:badInput',...
        'All inputs must be numbers!');
end
if(~isequalsizeorscalar(x,mu,sigma))
    error('seizmo:gaussian:badInput',...
        'All inputs must be equally sized or be scalars!');
end

% get pdf values
p=(1./(sigma*sqrt(2*pi))).*exp(-(x-mu).^2./(2*sigma.^2));

end
