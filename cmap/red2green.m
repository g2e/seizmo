function [map]=red2green(m)
%RED2GREEN    Red-White-Green colormap
%
%    Usage:    map=red2green(m)
%
%    Description:
%     MAP=RED2GREEN(M) returns a Mx3 matrix of RGB color values beginning
%     with red going to white in the middle and green at the end. RED2GREEN
%     by itself sets M to match the current figure's colormap size.  If no
%     figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(red2green)
%
%     % To go green to red:
%     colormap(flipud(red2green))
%
%    See also: BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN, DRYWET,
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
n=ceil(0.5*m);
if(mod(m,2))
    map=[ones(n,1) (0:n-1)'/(n-1) (0:n-1)'/(n-1); ...
        (m-n-1:-1:0)'./(m-n) ones(m-n,1) (m-n-1:-1:0)'./(m-n)];
else
    map=[ones(n,1) (0:n-1)'/n (0:n-1)'/n; ...
        (m-n-1:-1:0)'./(m-n) ones(m-n,1) (m-n-1:-1:0)'./(m-n)];
end

end