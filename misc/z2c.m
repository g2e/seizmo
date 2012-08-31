function [c]=z2c(z,cmap,clim)
%Z2C    Converts z-values to a color array of given color mapping & limits
%
%    Usage:    c=z2c(z,cmap)
%              c=z2c(z,cmap,clim)
%
%    Description:
%     C=Z2C(Z,CMAP) linearly maps Z values to colormap CMAP.  The lowest Z
%     value is mapped to the 1st color in CMAP, the maximum Z value is
%     mapped to last color in CMAP, and the remainder are given a color
%     based on their position within the range of Z and the number of
%     colors in CMAP.  Z must be a real valued vector.  CMAP must be a Mx3
%     array of RGB triplets where M is the number of colors in the mapping
%     spaced linearly between ZMIN & ZMAX.  C is a Nx3 array where N is the
%     number of elements in Z.  No color interpolation is done.
%
%     C=Z2C(Z,CMAP,CLIM) maps Z values to colormap CMAP using the map
%     limits CLIM rather than scaling to the limits of Z.  CLIM is expected
%     to be a 1x2 vector of [ZMIN ZMAX] over which the colormap changes.
%
%    Notes:
%     - NaN inputs are always set to [.25 .25 .25] which is a dark grey.
%
%    Examples:
%     % Plot some data with coloring based on their azimuth:
%     scatter(rand(1e3,1),rand(1e3,1),[],...
%             z2c(360*rand(1e3,1),hsv,[0 360]),...
%             'filled','markeredgecolor','k');
%
%    See also: NAME2RGB, COLORMAP, COLORMAPEDITOR

%     Version History:
%        Mar.  6, 2011 - initial version
%        Apr. 13, 2011 - fix bug when clim has no range
%        Aug. 30, 2012 - nan support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% default clim
if(nargin==2 || isempty(clim)); clim=[min(z) max(z)]; end

% fix nan color limits
if(all(isnan(clim))); clim=[min(z) max(z)]; end
if(all(isnan(clim))); clim=[0 1]; end

% fix equal color limits
if(diff(clim)==0); clim=clim(1)+eps(clim(1)).*[-1 1]; end

% check inputs
if(~isreal(z) || ~isvector(z))
    error('seizmo:z2c:badInput',...
        'Z must be a real-valued vector!');
elseif(ndims(cmap)~=2 || size(cmap,2)~=3 || any(cmap(:)<0 | cmap(:)>1))
    error('seizmo:z2c:badInput',...
        'CMAP must be a Nx3 matrix of valid rgb-triplets!');
elseif(~isreal(clim) || ~isvector(clim) || numel(clim)~=2)
    error('seizmo:z2c:badInput',...
        'CLIM must be a 1x2 vector of [ZMIN ZMAX]!');
end

% get color
c=.25*ones(numel(z),3);
nans=isnan(z);
nc=size(cmap,1);
idx=fix((z-clim(1))./(diff(clim)/nc));
idx(idx<1)=1;
idx(idx>nc)=nc;
c(~nans,:)=cmap(idx(~nans),:);

end
