function [corr]=geomsprcor(dd)
%GEOMSPRCOR    Geometrical spreading factor for amplitudes
%
%    Usage:    corr=geomsprcor(dd)
%
%    Description:
%     CORR=GEOMSPRCOR(DD) returns the geometrical spreading factor CORR for
%     each degree distance given in DD.  The formula for geometrical
%     spreading is rather simple:
%      CORR=(SIN(DD))^(-1/2).
%     The correction should be applied to amplitudes as follows:
%      Acorr=A/CORR.
%     This corrects all amplitudes to 90deg distance amplitude.
%
%    Notes:
%
%    Examples:
%     % Plot the spreading correction vs distance:
%     x=0:180;
%     y=geomsprcor(x);
%     figure; plot(x,y);
%     xlabel('degree distance')
%     ylabel('geometrical spreading factor')
%
%    See also: MANCOR, CRUCOR, ELLCOR

%     Version History:
%        Sep. 28, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2010 at 10:35 GMT

% todo

% check nargin
error(nargchk(1,1,nargin));

% check inputs
if(~isreal(dd) || any(dd<0))
    error('seizmo:geomsprcor:badInput',...
        'DD must be real-valued positive degree distance!');
end

% get spreading correction
corr=(sind(dd)).^(-1/2);

end
