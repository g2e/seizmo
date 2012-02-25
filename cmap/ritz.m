function [map]=ritz(m)
%RITZ    Colormap of Mike Ritzwoller's seismology group
%
%    Usage:    map=ritz(m)
%
%    Description:
%     MAP=RITZ(M) returns a Mx3 matrix of RGB color values going from black
%     to red, pink, orange, yellow, white, lightgreen, bluegray, cyan,
%     blue, & magenta.  RITZ by itself sets M to match the current figure's
%     colormap size.  If no figure exists, one is created.
%
%    Notes:
%     - Knockoff of Mike Ritzwoller's seismology group colormap
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(ritz)
%
%     % To use the colormap in reverse:
%     colormap(flipud(ritz))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, TOPO_LAND, FIRE,
%              NIGHTTIME, DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Feb. 17, 2010 - initial version
%        Mar. 22, 2010 - made bluegrey darker to balance better with orange
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 13:15 GMT

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x=(0:0.1:1)';
map=[interp1q(x,[0 0 0    1 2 3 4   6    6 8 10]'/12,(0:m-1)'/(m-1)) ...
     interp1q(x,[0 1 0.25 1 1 0 0.5 0.25 1 1 1 ]',(0:m-1)'/(m-1)) ...
     interp1q(x,[0 1 1    1 1 1 1   0.75 1 1 1 ]',(0:m-1)'/(m-1))];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
