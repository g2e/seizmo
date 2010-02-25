function [map]=ocean(m)
%OCEAN    Colormap for bathymetry
%
%    Usage:    map=ocean(m)
%
%    Description: MAP=OCEAN(M) returns a Mx3 matrix of RGB color values
%     going from black to seablue to white with strong hints of violet and
%     green along the way.  OCEAN by itself sets M to match the current
%     figure's colormap size.  If no figure exists, one is created.
%
%    Notes:
%     - by P. Wessel and F. Martinez, SOEST.
%     - The original GMT color palette extent was from -8000m to 0m.
%
%    Examples:
%     Set the current figure's colormap:
%      colormap(ocean)
%
%     To use the colormap in reverse:
%      colormap(flipud(ocean))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, DRYWET,
%              GEBCO, SEALAND, RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ

%     Version History:
%        Feb. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 17, 2010 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n=(0:m-1)'/(m-1);
x=(0:1/8:1)';
map=[interp1q(x,[0 0 0 0 0 86 172 211 250]'/255,n) ...
     interp1q(x,[0 5 10 80 150 197 245 250 255]'/255,n) ...
     interp1q(x,[0 25 50 125 200 184 168 211 255]'/255,n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
