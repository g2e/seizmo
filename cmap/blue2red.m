function [map]=blue2red(m)
%BLUE2RED    Blue-White-Red (aka Polar) colormap
%
%    Usage:    map=blue2red(m)
%
%    Description: MAP=BLUE2RED(M) returns a Mx3 matrix of RGB color values
%     beginning with blue going to white in the middle and red at the end.
%     BLUE2RED by itself sets M to match the current figure's colormap
%     size.  If no figure exists, one is created.
%
%    Notes:
%     - Matlab use to have a similar polar color map
%
%    Examples:
%     Set the current figure's colormap:
%      colormap(blue2red)
%
%     To go red to blue:
%      colormap(flipud(blue2red))
%
%    See also: RED2GREEN, GREEN2BLUE, SPLIT, SEIS, OCEAN, DRYWET,
%              GEBCO, SEALAND, RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ

%     Version History:
%        Feb. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 17, 2010 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n=ceil(0.5*m);
if(mod(m,2))
    map=[(0:n-1)'/(n-1) (0:n-1)'/(n-1) ones(n,1); ...
        ones(m-n,1) (m-n-1:-1:0)'./(m-n) (m-n-1:-1:0)'./(m-n)];
else
    map=[(0:n-1)'/n (0:n-1)'/n ones(n,1); ...
        ones(m-n,1) (m-n-1:-1:0)'./(m-n) (m-n-1:-1:0)'./(m-n)];
end

end
