function [x,t0]=swell_backproj(T,t)
%SWELL_BACKPROJ    Estimates the distance and origin time of a swell source
%
%    Usage:    [x,t0]=swell_backproj(T,t)
%
%    Description:
%     [x,t0]=SWELL_BACKPROJ(T,t) finds the best fitting distance and origin
%     time of a swell source given the swell info at a specific location
%     as periods T (in sec) and arrival times t (in sec).  The inputs must
%     contain at least 2 periods and arrival times in order to find a
%     distance and origin time of swell source.
%
%    Notes:
%     - Inputs should be equally sized or scalar.
%
%    Examples:
%     % Calculate the distance of a storm if 15s
%     % swell arrive 24hrs before the 10s swell:
%     x=swell_backproj([10 15],[0 -24*3600])
%
%     % Quick check:
%     swell_forward(10,x)-swell_forward(15,x)
%
%    See also: SWELL_FORWARD, SWELL_VELOCITY

%     Version History:
%        May  10, 2012 - initial version
%        May  16, 2012 - minor doc update, better checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  16, 2012 at 15:05 GMT

% todo:

error(nargchk(2,2,nargin));
if(~isscalar(T) && ~isscalar(t) && ~isequal(size(T),size(t)))
    error('seizmo:swell_forward:badInput',...
        'T & t must be equally sized or scalars!');
end
G=[4*pi/9.81./T(:) ones(numel(T),1)];
Gg=(G.'*G)\G.';
m=Gg*t(:);
x=m(1);
t0=m(2);

end
