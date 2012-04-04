function [cx,cy]=circle(r,steps)
%CIRCLE    Returns points on a circle in cartesian space
%
%    Usage:    [x,y]=circle(radius)
%              [x,y]=circle(radius,steps)
%
%    Description:
%     [X,Y]=CIRCLE(RADIUS) returns the cartesian coordinates of points on a
%     circle of radius RADIUS.  RADIUS must be a real-valued scalar.  The
%     points are spaced at 1 degree intervals.
%     
%     [X,Y]=CIRCLE(RADIUS,STEPS) sets the number of steps in the circle
%     given by X & Y.  The number of points in the circle is STEPS+1.
%     STEPS must be a positive integer >0.  The default is 360, which gives
%     a angular spacing of 1 degree.
%
%    Notes:
%
%    Examples:
%     % Draw a unit circle with radial lines every 30deg:
%     steps=360/30;
%     [x0,y0]=circle(0,steps);
%     [x1,y1]=circle(1,steps);
%     [x,y]=circle(1);
%     figure; plot([x0; x1],[y0; y1],'k',x,y,'k');
%     axis equal
%
%    See also: SPHERE, ELLIPSOID, CYLINDER

%     Version History:
%        May   4, 2010 - initial version
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 10:30 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default steps
if(nargin==1 || isempty(steps)); steps=360; end

% check input
if(~isreal(r) || ~isscalar(r))
    error('seizmo:circle:badInput',...
        'RADIUS must be a real valued scalar!');
elseif(~isreal(steps) || ~isscalar(steps) || fix(steps)~=steps || steps<1)
    error('seizmo:circle:badInput',...
        'STEPS must be an integer greater than 0!');
end

% get angles
ang=0:360/steps:360;

% get xy
cx=r*sind(ang);
cy=r*cosd(ang);

end
