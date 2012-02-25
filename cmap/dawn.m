function [map]=dawn(m)
%DAWN    Black-Blue-Violet-Pink-Orange-Red colormap
%
%    Usage:    map=dawn(m)
%
%    Description:
%     MAP=DAWN(M) returns a Mx3 matrix of RGB color values going through
%     black, blue, violet, pink, orange & red.  This colormap mimics the
%     array of colors in the gaining light before sunrise and has a strong
%     pastel scheme.  DAWN by itself sets M to match the current figure's
%     colormap size.  If no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(dawn)
%
%    See also: BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN, DRYWET,
%              NIGHTTIME, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ,
%              RED2GREEN, FIRE, DUSK, HSVSPIN, HSVCUSTOM

%     Version History:
%        July 25, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 01:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x=(0:1/6:1)';
map=[interp1q(x,[0 16 112 250 255 255 249]'./255,(0:m-1)'/(m-1)) ...
     interp1q(x,[0 20 109 191 178 161  88]'./255,(0:m-1)'/(m-1)) ...
     interp1q(x,[0 31 140 159 103  75  67]'./255,(0:m-1)'/(m-1))];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
