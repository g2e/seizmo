function [map]=sealand(m)
%SEALAND    Colormap for ocean and land
%
%    Usage:    map=sealand(m)
%
%    Description:
%     MAP=SEALAND(M) returns a Mx3 matrix of RGB color values going from
%     violet to cyan to green to peach to pink.  SEALAND by itself sets M
%     to match the current figure's colormap size.  If no figure exists,
%     one is created.
%
%    Notes:
%     - by W.H.F. Smith, NOAA
%     - The original GMT color palette extent was from -6000m to 3000m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(sealand)
%
%     % To use the colormap in reverse:
%     colormap(flipud(sealand))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, GMT_RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
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
x=[(0:12)'; (12:15)'; (15:18)']/18;
map=[interp1q(x,[255 240 225 210 195 180 165 150 135 120 105 90 75 ...
                 60 40 20 0 360 345 330 315]'/360,n) ...
     interp1q(x,[0.6*ones(1,13) 0.35 0.35 0.35 0.35 0.35 0.3 0.25 0.2]',n) ...
     ones(m,1)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
