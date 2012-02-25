function [map]=split(m)
%SPLIT    Lightblue-Blue-Black-Red-Lightred colormap
%
%    Usage:    map=split(m)
%
%    Description:
%     MAP=SPLIT(M) returns a Mx3 matrix of RGB color values going from
%     light blue to blue to black to red to light red.  SPLIT by itself
%     sets M to match the current figure's colormap size.  If no figure
%     exists, one is created.
%
%    Notes:
%     - matches GMT split colormap
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(split)
%
%     % To go in reverse:
%     colormap(flipud(split))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SEIS, GMT_OCEAN, DRYWET,
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
    map=[max((n-1:-1:0)'/(n-1)-0.5,0) max((n-1:-1:0)'/(n-1)-0.5,0) ...
        (n-1:-1:0)'/(n-1); (1:m-n)'/(m-n) ...
        max((1:m-n)'/(m-n)-0.5,0) max((1:m-n)'/(m-n)-0.5,0)];
else
    map=[max((n:-1:1)'/n-0.5,0) max((n:-1:1)'/n-0.5,0) ...
        (n:-1:1)'/n; (1:m-n)'/(m-n) ...
        max((1:m-n)'/(m-n)-0.5,0) max((1:m-n)'/(m-n)-0.5,0)];
end

end
