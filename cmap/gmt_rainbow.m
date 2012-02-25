function [map]=gmt_rainbow(m)
%GMT_RAINBOW    Standard rainbow colormap going from magenta to red
%
%    Usage:    map=gmt_rainbow(m)
%
%    Description:
%     MAP=GMT_RAINBOW(M) returns a Mx3 matrix of RGB color values going
%     from magenta to blue, cyan, green, yellow, orange, & red.
%     GMT_RAINBOW by itself sets M to match the current figure's colormap
%     size.  If no figure exists, one is created.
%
%    Notes:
%     - GMT_RAINBOW is not cyclic like HSV!
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(gmt_rainbow)
%
%     % To use the colormap in reverse:
%     colormap(flipud(gmt_rainbow))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ, FIRE, NIGHTTIME,
%              DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Feb. 17, 2010 - initial version
%        Feb. 22, 2012 - put a gmt_ to keep from conflicts
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
map=[interp1q([0; 1],[5/6; 0],(0:m-1)'/(m-1)) ones(m,1) ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
