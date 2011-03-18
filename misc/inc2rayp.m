function [rayp]=inc2rayp(inc,v,r)
%INC2RAYP    Returns the rayparameter for a given takeoff angle (from down)
%
%    Usage:    rayp=inc2rayp(inc,v,r)
%
%    Description:
%     RAYP=INC2RAYP(INC,V,R) returns the rayparameter RAYP for a ray in a
%     spherically symmetric Earth given that the ray has an inclination INC
%     at radius R and the velocity is V.  INC must be in degrees, R is in
%     km, V is in km/s.  RAYP is returned in s/deg.
%
%    Notes:
%
%    Examples:
%     % What is the ray parameter for a P-wave ray taking
%     % off at 15deg from a 600km deep earthquake:
%     model=prem('depth',600);
%     rayp=inc2rayp(15,model.vp,6371-600);
%
%    See also: RAYP2INC, RADPAT

%     Version History:
%        Mar.  8, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% get the rayparameter
rayp=pi/180*r.*sind(inc)./v;

end
