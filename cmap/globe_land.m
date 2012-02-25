function [map]=globe_land(m)
%GLOBE_LAND    Colormap for global topography
%
%    Usage:    map=globe_land(m)
%
%    Description:
%     MAP=GLOBE_LAND(M) returns a Mx3 matrix of RGB color values going from
%     green to peach, brown, light purple, & white. GLOBE_LAND by itself
%     sets M to match the current figure's colormap size.  If no figure
%     exists, one is created.
%
%    Notes:
%     - by Lester M. Anderson (CASP, UK)
%     - The original land portion of the GMT color palette extent was from
%       0m to 10000m.
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(globe_land)
%
%     % To use the colormap in reverse:
%     colormap(flipud(globe_land))
%
%    See also: RED2GREEN, BLUE2RED, GREEN2BLUE, SPLIT, SEIS, GMT_OCEAN,
%              DRYWET, GEBCO, SEALAND, GMT_RAINBOW, RELIEF, GLOBE_SEA,
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
x=[0; 0.01; 0.02; (1:20)'/20];
map=[interp1q(x,[51 51 187 255 243 230 217 168 164 162 159 156 ...
                 153 162 178 183 194 204 229 242 255 255 255]'/255,n) ...
     interp1q(x,[102 204 228 220 202 184 166 154 144 134 123 113 ...
                 102 89 118 147 176 204 229 242 255 255 255]'/255,n) ...
     interp1q(x,[0 102 146 185 137 88 39 31 25 19 13 7 0 89 118 ...
                 147 176 204 229 242 255 255 255]'/255,n)];

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

end
