function [x]=triangletf(t,t0,hwidth,amp)
%TRIANGLETF    Returns a triangle time function
%
%    Usage:    x=triangletf(t,t0,hwidth)
%              x=triangletf(t,t0,hwidth,amp)
%
%    Description:
%     X=TRIANGLETF(T,T0,HWIDTH) creates a triangle time function centered
%     at T0 that is 0 at times at/outside T0-HWIDTH and T0+HWIDTH.  The
%     peak value at T0 is 1/HWIDTH so that the area under the curve is 1.
%
%     X=TRIANGLETF(T,T0,HWIDTH,AMP) sets the peak amplitude of the triangle
%     (at time T0) to AMP.  The default value of AMP is 1/HWIDTH.
%
%    Notes:
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
%    See also: GAUSSIANTF, TAPERFUN, TRIANG

%     Version History:
%        Oct. 17, 2009 - initial version
%        Feb.  1, 2011 - default amp sets area to 1, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  1, 2011 at 21:15 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% defaults
if(nargin<4); amp=[]; end

% check inputs
if(~isnumeric(t))
    error('seizmo:gaussiantf:badInput','T must be a numeric array!');
elseif(~isscalar(t0) || ~isnumeric(t0))
    error('seizmo:gaussiantf:badInput','TO must be a numeric scalar!');
elseif(~isscalar(hwidth) || ~isnumeric(hwidth))
    error('seizmo:gaussiantf:badInput','HWIDTH must be a numeric scalar!');
elseif(~isempty(amp) && (~isscalar(amp) || ~isnumeric(amp)))
    error('seizmo:gaussiantf:badInput','AMP must be a numeric scalar!');
end

% get triangle values
if(isempty(amp)); amp=1/hwidth; end
x=amp.*(1-abs((t-t0)./hwidth));
x(x<0)=0;

end
