function [n,e,u]=strikedip2neu(strike,dip)
%STRIKEDIP2NEU    Returns the normal vector to the plane given strike & dip
%
%    Usage:    neu=strikedip2neu(sd)
%              neu=strikedip2neu(strike,dip)
%              [n,e,u]=strikedip2neu(...)
%
%    Description:
%     NEU=STRIKEDIP2NEU(SD) converts the strike & dip of a fault plane into
%     the vector normal to that fault as [North East Up] in NEU.  SD is
%     expected as [strike dip] in degrees where strike is clockwise
%     from North and dip is positive downward from the horizontal.  Note
%     that the strike must be such that when you look along the direction
%     of the strike the fault dips to your right.  SD may be a Nx2 array to
%     allow for conversion of multiple fault planes (NEU is Nx3).
%
%     NEU=STRIKEDIP2NEU(STRIKE,DIP) allows strike & dip to be given
%     separately.
%
%     [N,E,U]=STRIKEDIP2NEU(...) returns North, East & Up separately.
%
%    Notes:
%     - The returned normal vector always points towards the hanging wall.
%
%    Examples:
%     % What is the plunge & azimuth of the normal vector
%     % to a fault striking North and dipping at 30deg?
%     neu2vpa(strikedip2neu(0,30))
%
%     % "Fix" the upwards normal vector:
%     neu2vpa(-1*strikedip2neu(0,30))
%
%    See also: NEU2STRIKEDIP, AUXPLANE

%     Version History:
%        Mar. 14, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 14, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% one or both inputs
switch nargin
    case 1
        if(size(strike,2)~=2 || ndims(strike)>2)
            error('seizmo:strikedip2neu:badInput',...
                'SD must be a Nx2 array as [STRIKE DIP] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:strikedip2neu:badInput',...
                'SD must be a real-valued Nx2 array!');
        end
        [strike,dip]=deal(strike(:,1),strike(:,2));
    case 2
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip))
            error('seizmo:strikedip2neu:badInput',...
                'STRIKE/DIP must be real-valued arrays!');
        end
        [strike,dip]=expandscalars(strike,dip);
end

% calculate north/east/up
[n,e,u]=deal(-sind(dip).*sind(strike),sind(dip).*cosd(strike),cosd(dip));

% combine if only one output
if(nargout<=1); n=[n(:) e(:) u(:)]; end

end
