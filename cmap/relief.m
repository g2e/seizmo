function [map]=relief(m)
%RELIEF    Colormap for global topography
%
%    Usage:    map=relief(m)
%
%    Description:
%     MAP=RELIEF(M) returns a Mx3 matrix of RGB color values going from
%     green quickly to brown and then to white.  RELIEF by itself sets M to
%     match the current figure's colormap size.  If no figure exists, one
%     is created.
%
%    Notes:
%     - by P. Wessel and F. Martinez, SOEST
%     - The original land portion of the GMT color palette extent was from
%       0m to 8000m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(relief)
%
%     % To use the colormap in reverse:
%     colormap(flipud(relief))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, GLOBE_SEA, GLOBE_LAND,
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
x=[0; 1/16; (1:8)'/8];
map=[interp1q(x,[70 120 146 198 250 250 252 252 253 255]'/255,n) ...
     interp1q(x,[120 100 126 178 230 234 238 243 249 255]'/255,n) ...
     interp1q(x,[50 50 60 80 100 126 152 177 216 255]'/255,n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
