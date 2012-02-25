function [map]=fire(m)
%FIRE    Black-Blue-Violet-Red-Yellow-White colormap
%
%    Usage:    map=fire(m)
%
%    Description:
%     MAP=FIRE(M) returns a Mx3 matrix of RGB color values going from black
%     to blue to violet to red to yellow to white.  This colormap can make
%     an eerie night fire.  It is used primarily for spectrogram plots.
%     FIRE by itself sets M to match the current figure's colormap size.
%     If no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(fire)
%
%    See also: BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN, DRYWET,
%              NIGHTTIME, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
%              GLOBE_LAND, SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ,
%              RED2GREEN, DUSK, DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        Apr. 26, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 01:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x=(0:0.2:1)';
map=[interp1q(x,[0 0 1 1 1 1]',(0:m-1)'/(m-1)) ...
     interp1q(x,[0 0 0 0 1 1]',(0:m-1)'/(m-1)) ...
     interp1q(x,[0 1 1 0 0 1]',(0:m-1)'/(m-1))];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
