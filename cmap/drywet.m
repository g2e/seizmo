function [map]=drywet(m)
%DRYWET    Colormap for hydrology
%
%    Usage:    map=drywet(m)
%
%    Description:
%     MAP=DRYWET(M) returns a Mx3 matrix of RGB color values going from
%     brown to light blue to a deep blue.  DRYWET by itself sets M to match
%     the current figure's colormap size.  If no figure exists, one is
%     created.
%
%    Notes:
%     - by Ed Maurer, U Washington
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(drywet)
%
%     % To use the colormap in reverse:
%     colormap(flipud(drywet))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ,
%              FIRE, NIGHTTIME, DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Feb. 17, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n=(0:m-1)'/(m-1);
x=(0:1/6:1)';
map=[interp1q(x,[134 238 180 50 12 38 8]'/255,n) ...
     interp1q(x,[97 199 238 238 120 1 51]'/255,n) ...
     interp1q(x,[42 100 135 235 238 183 113]'/255,n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
