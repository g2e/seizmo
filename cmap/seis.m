function [map]=seis(m)
%SEIS    Colormap for seismology
%
%    Usage:    map=seis(m)
%
%    Description:
%     MAP=SEIS(M) returns a Mx3 matrix of RGB color values going from red
%     to orange to yellow to green to blue.  SEIS by itself sets M to match
%     the current figure's colormap size.  If no figure exists, one is
%     created.
%
%    Notes:
%     -  by Susan van der Lee
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(seis)
%
%     % To use the colormap in reverse:
%     colormap(flipud(seis))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, GMT_OCEAN, DRYWET,
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
n=(0:m-1)'/(m-1);
x=(0:1/9:1)';
map=[interp1q(x,[2/3 1 1 1 1 1 6/17 0 0 0]',n) ...
     interp1q(x,[0 0 1/3 2/3 1 1 1 16/17 16/51 0]',n) ...
     interp1q(x,[0 0 0 0 0 0 2/17 22/51 1 41/51]',n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
