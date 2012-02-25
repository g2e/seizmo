function [map]=sealand_land(m)
%SEALAND_LAND    Land portion of the sealand colormap
%
%    Usage:    map=sealand_land(m)
%
%    Description:
%     MAP=SEALAND_LAND(M) returns a Mx3 matrix of RGB color values going
%     from light yellow to light red to pink.  SEALAND_LAND by itself sets
%     M to match the current figure's colormap size.  If no figure exists,
%     one is created.
%
%    Notes:
%     - by W.H.F. Smith, NOAA
%     - The original land portion of the GMT color palette extent was from
%       0m to 3000m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(sealand_land)
%
%     % To use the colormap in reverse:
%     colormap(flipud(sealand_land))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, TOPO_LAND, RITZ, FIRE, NIGHTTIME,
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
x=[(0:3)'; (3:6)']/6;
map=[interp1q(x,[60 40 20 0 360 345 330 315]'/360,n) ...
     interp1q(x,[0.35 0.35 0.35 0.35 0.35 0.3 0.25 0.2]',n) ...
     ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
