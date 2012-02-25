function map=hsvcustom(hsv,m)
%HSVCUSTOM    HSV colormap with start/end at color specified
%
%    Usage:    map=hsvcustom(hsv,m)
%
%    Description:
%     MAP=HSVCUSTOM(HSV,M) returns a Mx3 matrix of RGB color values
%     corresponding to varying hue in the hue-saturation-value color
%     system.  This function differs from HSV in that it will begin with
%     the hue specified by HSV (0 to 1) rather than always starting/ending
%     with red).  The saturation and value are not changed across the color
%     map.  Calling HSVCUSTOM with no inputs is equivalent to HSV.  Not
%     specifying M sets it to match the current figure's colormap size.  If
%     no figure exists, one is created.
%
%    Notes:
%
%    Examples:
%     % Set the current figure's colormap:
%     colormap(hsvcustom([0.3 0.5 0.5]))
%
%     % To use the colormap in reverse:
%     colormap(flipud(hsvcustom([0.3 0.5 0.5])))
%
%    See also: HSV, HSV2RGB, RGB2HSV, RED2GREEN, BLUE2RED, GREEN2BLUE,
%              SPLIT, SEIS, GMT_OCEAN, DRYWET, GEBCO, SEALAND, GMT_RAINBOW,
%              RELIEF, GLOBE_SEA, GLOBE_LAND, SEALAND_SEA, SEALAND_LAND,
%              TOPO_LAND, FIRE, HSVSPIN, NIGHTTIME, DUSK, DAWN

%     Version History:
%        May  11, 2010 - initial version
%        Feb. 22, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2012 at 00:45 GMT

if nargin < 1, hsv=[0 1 1]; end
if nargin < 2, m = size(get(gcf,'colormap'),1); end
h = mod(hsv(1)+(0:m-1)'/max(m,1),1);
if isempty(h)
  map = [];
else
  map = hsv2rgb([h hsv(2)*ones(m,1) hsv(3)*ones(m,1)]);
end
map(map>1)=1;
map(map<0)=0;

end
