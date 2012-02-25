function [map]=green2blue(m)
%GREEN2BLUE    Green-White-Blue colormap
%
%    Usage:    map=green2blue(m)
%
%    Description:
%     MAP=GREEN2BLUE(M) returns a Mx3 matrix of RGB color values beginning
%     with green going to white in the middle and blue at the end.
%     GREEN2BLUE by itself sets M to match the current figure's colormap
%     size.  If no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(green2blue)
%
%     % To go blue to green:
%     colormap(flipud(green2blue))
%
%    See also: RED2GREEN, BLUE2RED, SPLIT, SEIS, GMT_OCEAN, DRYWET,
%              GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ, FIRE, NIGHTTIME,
%              DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Feb. 17, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n=ceil(0.5*m);
if(mod(m,2))
    map=[(0:n-1)'/(n-1) ones(n,1) (0:n-1)'/(n-1); ...
        (m-n-1:-1:0)'./(m-n) (m-n-1:-1:0)'./(m-n) ones(m-n,1)];
else
    map=[(0:n-1)'/n ones(n,1) (0:n-1)'/n; ...
        (m-n-1:-1:0)'./(m-n) (m-n-1:-1:0)'./(m-n) ones(m-n,1)];
end

end
