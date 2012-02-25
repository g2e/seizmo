function [map]=dusk(m)
%DUSK    Black-Blue-Violet-Pink-Tan-White colormap
%
%    Usage:    map=dusk(m)
%
%    Description:
%     MAP=DUSK(M) returns a Mx3 matrix of RGB color values going through
%     black, blue, violet, pink, tan & white.  This colormap mimics the
%     array of colors in the fading light after sunset and has a strong
%     pastel scheme.  DUSK by itself sets M to match the current figure's
%     colormap size.  If no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(dusk)
%
%    See also: BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN, DRYWET,
%              NIGHTTIME, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ,
%              RED2GREEN, FIRE, DAWN

%     Version History:
%        July 25, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 01:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x=(0:1/6:1)';
map=[interp1q(x,[0 16 134 255 243 255 255]'./255,(0:m-1)'/(m-1)) ...
     interp1q(x,[0 20 116 151 168 201 255]'./255,(0:m-1)'/(m-1)) ...
     interp1q(x,[0 31 157 158 147 152 255]'./255,(0:m-1)'/(m-1))];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
