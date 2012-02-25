function map=hsvspin(m)
%HSVSPIN    HSV colormap with random circular shift
%
%    Usage:    map=hsvspin(m)
%
%    Description:
%     MAP=HSVSPIN(M) returns a Mx3 matrix of RGB color values corresponding
%     to varying hue in the hue-saturation-value color system.  This
%     function differs from HSV in that it will begin with a random fully
%     saturated color (rather than always starting/ending with red).
%     HSVSPIN by itself sets M to match the current figure's colormap
%     size. If no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(hsvspin)
%
%     % To use the colormap in reverse:
%     colormap(flipud(hsvspin))
%
%    See also: HSV, HSV2RGB, RGB2HSV, RED2GREEN, BLUE2RED, GREEN2BLUE,
%              SPLIT, SEIS, GMT_OCEAN, DRYWET, GEBCO, SEALAND, GMT_RAINBOW,
%              RELIEF, GLOBE_SEA, GLOBE_LAND, SEALAND_SEA, SEALAND_LAND,
%              TOPO_LAND, FIRE, NIGHTTIME, DUSK, DAWN, HSVSPIN

%     Version History:
%        May   7, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:45 GMT

if nargin < 1, m = size(get(gcf,'colormap'),1); end
h = mod(rand+(0:m-1)'/max(m,1),1);
if isempty(h)
  map = [];
else
  map = hsv2rgb([h ones(m,2)]);
end

end
