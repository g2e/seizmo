function [n,e,u]=strikedip2norm(strike,dip,varargin)
%STRIKEDIP2NORM    Returns the normal vector to a fault plane
%
%    Usage:    neu=strikedip2norm(sd)
%              neu=strikedip2norm(strike,dip)
%              [n,e,u]=strikedip2norm(...)
%
%    Description:
%     NEU=STRIKEDIP2NORM(SD) converts the strike & dip of a fault plane to
%     the vector normal to the fault as [North East Up] in NEU.  SD is
%     expected as [strike dip] in degrees where strike is clockwise
%     from North and dip is positive downward from the horizontal.  Note
%     that the strike must be such that when you look along the direction
%     of the strike the fault dips to your right.  SD may be a Nx2 array to
%     allow for conversion of multiple fault planes (NEU is Nx3).
%
%     NEU=STRIKEDIP2NORM(STRIKE,DIP) allows strike & dip to be given
%     separately.
%
%     [N,E,U]=STRIKEDIP2NORM(...) returns North, East & Up separately.
%
%    Notes:
%     - The returned normal vector always points towards the hanging wall.
%       This means the normal vector always points to the right when
%       looking along the strike (except for a no dip horizontal fault).
%
%    Examples:
%     % What is the plunge & azimuth of the normal vector
%     % to a fault striking North and dipping at 30deg?
%     neu2vpa(strikedip2norm(0,30))
%
%     % "Fix" the upwards normal vector so it points to the foot wall:
%     neu2vpa(-1*strikedip2norm(0,30))
%
%    See also: NORM2STRIKEDIP, SDR2NULL, SDR2SLIP, AUXPLANE, SDR2TPB,
%              TPB2SDR, NORMSLIP2SDR, NODALLINES

%     Version History:
%        Mar. 14, 2013 - initial version
%        Mar. 15, 2013 - rename for clarity
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% one or both inputs
switch nargin
    case 1
        % ignore the rake term if present
        if(~any(size(strike,2)==[2 3]) || ndims(strike)>2)
            error('seizmo:strikedip2norm:badInput',...
                'SD must be a Nx2 array as [STRIKE DIP] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:strikedip2norm:badInput',...
                'SD must be a real-valued Nx2 array!');
        end
        [strike,dip]=deal(strike(:,1),strike(:,2));
    case {2 3}
        % ignore the rake term if present
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip))
            error('seizmo:strikedip2norm:badInput',...
                'STRIKE/DIP must be real-valued arrays!');
        end
        [strike,dip]=expandscalars(strike,dip);
end

% calculate north/east/up
[n,e,u]=deal(-sind(dip).*sind(strike),sind(dip).*cosd(strike),cosd(dip));

% combine if only one output
if(nargout<=1); n=[n(:) e(:) u(:)]; end

end
