function [map]=globe_sea(m)
%GLOBE_SEA    Colormap for global bathymetry
%
%    Usage:    map=globe_sea(m)
%
%    Description:
%     MAP=GLOBE_SEA(M) returns a Mx3 matrix of RGB color values going from
%     violet to light blue to white.  GLOBE_SEA by itself sets M to match
%     the current figure's colormap size.  If no figure exists, one is
%     created.
%
%    Notes:
%     - by Lester M. Anderson (CASP, UK)
%     - The original sea portion of the GMT color palette extent was from
%       -10000m to 0m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(globe_sea)
%
%     % To use the colormap in reverse:
%     colormap(flipud(globe_sea))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_LAND,
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
x=[(0:19)'/20; 49/50; 1];
map=[interp1q(x,[153 153 153 136 119 102 85 68 51 34 17 ...
                 0 27 54 81 108 134 161 188 215 241 241]'/255,n) ...
     interp1q(x,[0 0 0 17 34 51 68 85 102 119 136 153 ...
                 164 175 186 197 208 219 230 241 252 252]'/255,n) ...
     ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
