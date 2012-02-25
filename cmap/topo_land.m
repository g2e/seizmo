function [map]=topo_land(m)
%TOPO_LAND    Colormap used for topography in tectonic maps
%
%    Usage:    map=topo_land(m)
%
%    Description:
%     MAP=TOPO_LAND(M) returns a Mx3 matrix of RGB color values going from
%     light blue quickly to light green and tan and eventually to white.
%     TOPO_LAND by itself sets M to match the current figure's colormap
%     size.  If no figure exists, one is created.
%
%    Notes:
%     - by D. Sandwell, Scripps
%     - The original land portion of the GMT color palette extent was from
%       0m to 7000m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(topo_land)
%
%     % To use the colormap in reverse:
%     colormap(flipud(topo_land))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, RITZ, FIRE,
%              NIGHTTIME, DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Feb. 17, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n=(0:m-1)'/(m-1);
x=[0; (1:3)'/35; 1/7; 3/14; 0.5; 1];
map=[interp1q(x,[195 160 125 99 75 50 25 0]'/360,n) ...
     interp1q(x,[0.35 0.4 0.45 0.45 0.45 0.35 0.1 0]',n) ...
     interp1q(x,[0.7 0.7 0.7 0.8 0.8 0.9 1 1]',n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
