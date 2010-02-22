function [map]=rainbow(m)
%RAINBOW    Standard rainbow color map going from magenta to red
%
%    Usage:    map=rainbow(m)
%
%    Description: MAP=RAINBOW(M) returns a Mx3 matrix of RGB color values
%     going from magenta to blue, cyan, green, yellow, orange, & red.
%     RAINBOW by itself sets M to match the current figure's colormap size.
%     If no figure exists, one is created.
%
%    Notes:
%     - RAINBOW is not cyclic like HSV!
%
%    Examples:
%     Set the current figure's colormap:
%      colormap(rainbow)
%
%     To use the colormap in reverse:
%      colormap(flipud(rainbow))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, OCEAN, DRYWET,
%              GEBCO, SEALAND, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ

%     Version History:
%        Feb. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 17, 2010 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
map=[interp1q([0; 1],[5/6; 0],(0:m-1)'/(m-1)) ones(m,1) ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
