function [x]=gaussiantf(t,t0,hwidth,hw2std,amp)
%GAUSSIANTF    Returns a gaussian time function
%
%    Usage:    x=gaussiantf(t,t0,hwidth)
%              x=gaussiantf(t,t0,hwidth,hw2std)
%              x=gaussiantf(t,t0,hwidth,hw2std,amp)
%
%    Description:
%     X=GAUSSIANTF(T,T0,HWIDTH) samples a gaussian curve centered at time
%     TO at the times in array T.  The gaussian curve is further defined to
%     falloff to a value of 1/e at a time of HWIDTH from TO.  Both TO and
%     HWIDTH must be numeric scalars.  T must be a numeric array and X is
%     an equal sized array containing the associated gaussian values.
%
%     X=GAUSSIANTF(T,T0,HWIDTH,HW2STD) scales the half width into standard
%     deviations.  The default value is exp(0.5)*sqrt(2), which creates a
%     gaussian that drops off similarly to a a triangular source time
%     function with the same halfwidth.  See the Notes section for the
%     gaussian formula.
%
%     X=GAUSSIANTF(T,T0,HWIDTH,HW2STD,AMP) sets the peak amplitude of the
%     gaussian curve at T0.  The default value of AMP is set so the area
%     under the gaussian curve sums to 1.
%
%    Notes:
%     - Formula for returned gaussian values:
%
%                      1     /           T-T0  \  2
%                   - ___ * | HW2STD * ________ |
%                      2     \          HWIDTH /
%        X = AMP * e
%
%       where HWIDTH/HW2STD = SIGMA (the width of one standard deviation)
%       and AMP = HW2STD/(HWIDTH*sqrt(2*PI)) which scales the gaussian so
%       the area under the curve is 1.
%
%    Examples:
%     % Compare gaussian and triangle source time functions
%     % for convolving with synthetic seismic data:
%     t=(-15:0.1:15)';
%     fh=figure; ax=axes('parent',fh);
%     plot(ax,t,[gaussiantf(t,0,10) triangletf(t,0,10)],'linewidth',2);
%     legend(ax,{'Gaussian' 'Triangular'});
%     xlabel(ax,'Time (s)');
%     grid(ax,'on');
%
%    See also: TRIANGLETF, TAPERFUN, GAUSSWIN

%     Version History:
%        Oct. 11, 2009 - initial version
%        Oct. 17, 2009 - new example
%        Feb.  1, 2011 - greatly improved docs, renamed variable that
%                        scaled hw to sigma, default amp sets area to 1,
%                        default hw2std set to mimic triangletf
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  1, 2011 at 21:30 GMT

% todo:

% check nargin
error(nargchk(3,5,nargin));

% defaults
if(nargin<4 || isempty(hw2std)); hw2std=exp(0.5)*sqrt(2); end
if(nargin<5); amp=[]; end

% check inputs
if(~isnumeric(t))
    error('seizmo:gaussiantf:badInput','T must be a numeric array!');
elseif(~isscalar(t0) || ~isnumeric(t0))
    error('seizmo:gaussiantf:badInput','TO must be a numeric scalar!');
elseif(~isscalar(hwidth) || ~isnumeric(hwidth))
    error('seizmo:gaussiantf:badInput','HWIDTH must be a numeric scalar!');
elseif(~isscalar(hw2std) || ~isnumeric(hw2std))
    error('seizmo:gaussiantf:badInput','HW2STD must be a numeric scalar!');
elseif(~isempty(amp) && (~isscalar(amp) || ~isnumeric(amp)))
    error('seizmo:gaussiantf:badInput','AMP must be a numeric scalar!');
end

% get gaussian values
if(isempty(amp)); amp=hw2std/(hwidth*sqrt(2*pi)); end
x=amp.*exp(-1./2.*(hw2std.*(t-t0)./hwidth).^2);

end
