function [map]=sealand_sea(m)
%SEALAND_SEA    Sea portion of the sealand colormap
%
%    Usage:    map=sealand_sea(m)
%
%    Description:
%     MAP=SEALAND_SEA(M) returns a Mx3 matrix of RGB color values going
%     from violet to cyan to green to yellow.  SEALAND_SEA by itself sets M
%     to match the current figure's colormap size.  If no figure exists,
%     one is created.
%
%    Notes:
%     - by W.H.F. Smith, NOAA
%     - The original sea portion of the GMT color palette extent was from
%       -6000m to 0m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(sealand_sea)
%
%     % To use the colormap in reverse:
%     colormap(flipud(sealand_sea))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_LAND, TOPO_LAND, RITZ, FIRE, NIGHTTIME,
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
x=(0:12)'/12;
map=[interp1q(x,(17:-1:5)'/24,n) 0.6*ones(m,1) ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
