function [map]=nighttime(m)
%NIGHTTIME    For NASA nighttime lights maps
%
%    Usage:    map=nighttime(m)
%
%    Description:
%     MAP=NIGHTTIME(M) returns a Mx3 matrix of RGB color values (even the
%     map is created in HSV space) going from a dark to light blue, then
%     sharply to tan grading to a light pink.  This colormap is useful when
%     plotting NASA nighttime datasets.  NIGHTTIME by itself sets M to
%     match the current figure's colormap size.  If no figure exists, one
%     is created.
%
%    Notes:
%     - originally written by Andreas Trawoeger to match NASA's colormap
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(nighttime)
%
%    See also: BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN, DRYWET, FIRE,
%              GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA, GLOBE_LAND,
%              SEALAND_SEA, SEALAND_LAND, TOPO_LAND, RITZ, RED2GREEN, DUSK,
%              DAWN, HSVSPIN, HSVCUSTOM

%     Version History:
%        July 25, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 01:40 GMT

% todo:

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x=[0 0.5 0.5+eps 1]';
map=[interp1q(x,[260  195   65 0.0]'./360,(0:m-1)'/(m-1)) ...
     interp1q(x,[1.0 0.55 0.55 0.1]',(0:m-1)'/(m-1)) ...
     interp1q(x,[0.1 0.55 0.55 1.0]',(0:m-1)'/(m-1))];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% convert to rgb
map=hsv2rgb(map);

end
