function [map]=gebco(m)
%GEBCO    Colormap approximating GEBCO bathymetry charts
%
%    Usage:    map=gebco(m)
%
%    Description:
%     MAP=GEBCO(M) returns a Mx3 matrix of RGB color values going from cyan
%     to light greens to white.  GEBCO by itself sets M to match the
%     current figure's colormap size.  If no figure exists, one is created.
%
%    Notes:
%     - by Andrew Goodwillie, Scripps
%     - The original GMT color palette extent was from -7000m to 0m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(gebco)
%
%     % To use the colormap in reverse:
%     colormap(flipud(gebco))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
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
n=(0:m-1)'/(m-1);
x=(0:6)'/7;
x=[x; 13/14; 34/35; 1];

map=[interp1q(x,[0 35 90 140 165 195 210 230 235 235]'/255,n) ...
     interp1q(x,[16/17 1 1 1 1 1 1 1 1 1]',n) ...
     interp1q(x,[255 255 255 230 215 215 215 240 255 255]'/255,n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
