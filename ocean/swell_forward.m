function [t]=swell_forward(T,x)
%SWELL_FORWARD    Predicts the time of ocean swell given period & distance
%
%    Usage:    t=swell_forward(T,x)
%
%    Description:
%     t=SWELL_FORWARD(T,x) predicts the travel time of swells through deep
%     ocean given the period T in seconds and distance traveled x in
%     meters.  The relationship is:
%
%           4*PI*x
%      t = ________
%            g*T
%
%     where g is the gravitational acceleration 9.81m/s^2.  The output is
%     in seconds.
%
%    Notes:
%     - Inputs should be equally sized or scalar.
%
%    Examples:
%     % Find the arrival time of 15s swell from a short storm 5000km away
%     % occuring at noon on March 15, 2012:
%     t0=datenum([2012 3 15 12 0 0]);
%     t=swell_forward(15,5000);
%     t=datevec(t/86400+t0);
%     datestr(t)
%
%    See also: SWELL_BACKPROJ, SWELL_VELOCITY

%     Version History:
%        May  10, 2012 - initial version
%        May  16, 2012 - minor doc update, better checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  16, 2012 at 15:05 GMT

% todo:

error(nargchk(2,2,nargin));
if(~isscalar(T) && ~isscalar(x) && ~isequal(size(T),size(x)))
    error('seizmo:swell_forward:badInput',...
        'T & x must be equally sized or scalars!');
end
t=4.*pi.*x./9.81./T;

end
