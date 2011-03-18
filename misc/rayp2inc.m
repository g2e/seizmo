function [inc]=rayp2inc(rayp,v,r)
%RAYP2INC    Returns the takeoff angle (from down) for a given rayparameter
%
%    Usage:    inc=rayp2inc(rayp,v,r)
%
%    Description:
%     INC=RAYP2INC(RAYP,V,R) calculates the inclination of a ray path in a
%     spherically symmetric Earth given the ray parameter RAYP and velocity
%     V at a specified radius R.  RAYP is in s/deg, V is in km/s, and
%     radius is in km.  INC is returned in degrees.
%
%    Notes:
%
%    Examples:
%     % Find the takeoff angle for a core-diffracted
%     % p-wave from a 600km deep earthquake:
%     model=prem('depth',600);
%     pdiff=tauptime('mod','prem','ph','Pdiff','deg',120);
%     takeoff=rayp2inc(pdiff.rayparameter,model.vp,6371-600);
%
%    See also: INC2RAYP, RADPAT

%     Version History:
%        Mar.  8, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% get the takeoff angle
inc=asind(180/pi*rayp.*v./r);

end
