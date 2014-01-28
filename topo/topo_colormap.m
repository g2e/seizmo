function [varargout]=topo_colormap(varargin)
%TOPO_COLORMAP    Creates a colormap for topography data
%
%    Usage:    topo_colormap(z)
%              topo_colormap(z,m)
%              topo_colormap(z,m,sea_cmap,land_cmap)
%              topo_colormap(ax,...)
%              [cmap,z0]=topo_colormap(...)
%
%    Description:
%     TOPO_COLORMAP(Z) creates and applies a colormap for the topography
%     data Z.  The colormap is a combination of two colormaps appropriately
%     proportioned to the elevation ranges of Z so that the interface
%     between the colormaps is at sea level.  The two default colormaps are
%     provided by the functions GLOBE_SEA & GLOBE_LAND.  With no output
%     arguments, the figure corresponding to the current axes is set to use
%     the created colormap and the color axis is adjusted to the limits of
%     Z (unless Z does not span 0 - see the Notes section below).
%
%     TOPO_COLORMAP(Z,M) uses a colormap with M colors.  The default is
%     that of the current figure's colormap size (the default is 64).
%
%     TOPO_COLORMAP(Z,M,SEA_CMAP,LAND_CMAP) uses the specified colormaps to
%     create the topography colormap.  SEA_CMAP & LAND_CMAP should be Nx3
%     arrays of rgb colors.  The default sea and land colormaps are sampled
%     at M colors, combined and resampled to M colors total.
%
%     TOPO_COLORMAP(AX,...) sets the colormap of the figure corresponding
%     to the axes AX.
%
%     [CMAP,Z0]=TOPO_COLORMAP(...) outputs the colormap as CMAP and the
%     z-values of the colormap as Z0.
%
%    Notes:
%     - The output colormap is always adjusted to contain 0.  This means
%       that Z0 will always span [min([Z(:); 0]) max([Z(:); 0])].  This
%       makes the absolute values of Z more apparent while reducing the
%       highlighting of the range of Z values.
%
%    Examples:
%     % Force the colormap to be from -7km to 5km with 256 colors:
%     cmap=topo_colormap([-7000 5000],256);
%
%    See also: TOPO_REGION, TOPO_POINTS

%     Version History:
%        Feb. 19, 2010 - initial version
%        Feb. 22, 2010 - output all z-values rather than limits
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  2, 2012 - minor doc update
%        Jan. 28, 2014 - added axes specification option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% parse axes
[ax,varargin]=axparse(varargin{:});
nargs=numel(varargin);

% check nargin
error(nargchk(1,4,nargs));

% assign variables
if(numel(varargin)<4); varargin(numel(varargin)+1:4)={[]}; end
[z,m,sea,land]=deal(varargin{:});

% default axes
if(isempty(ax)); ax=gca; end

% check inputs
if(isempty(m)); m=size(get(get(ax,'parent'),'colormap'),1); end
if(isempty(sea)); sea=globe_sea(m); end
if(isempty(land)); land=globe_land(m); end
if(isempty(z))
    error('seizmo:topo_colormap:badInput',...
        'Z must not be empty!');
elseif(~isnumeric(z))
    error('seizmo:topo_colormap:badInput',...
        'Z must be numeric!');
elseif(~isreal(z))
    warning('seizmo:topo_colormap:complexInput',...
        ['Z is complex.  Only the real-portion of Z\n' ...
        'will be used for the color map range!']);
elseif(~isscalar(m) || ~isreal(m) || m~=fix(m) || m<=0)
    error('seizmo:topo_colormap:badInput',...
        'M must be a positive scalar integer!');
elseif(ndims(sea)>2 || size(sea,2)~=3)
    error('seizmo:topo_colormap:badInput',...
        'SEA_CMAP must be a Nx3 matrix!');
elseif(~isreal(sea) || any(sea(:)>1 | sea(:)<0))
    error('seizmo:topo_colormap:badInput',...
        'SEA_CMAP must have real-valued elements in the range [0 1]!');
elseif(ndims(land)>2 || size(land,2)~=3)
    error('seizmo:topo_colormap:badInput',...
        'LAND_CMAP must be a Nx3 matrix!');
elseif(~isreal(land) || any(land(:)>1 | land(:)<0))
    error('seizmo:topo_colormap:badInput',...
        'LAND_CMAP must have real-valued elements in the range [0 1]!');
end

% get range of z
zmin=min([real(z(:)); 0]);
zmax=max([real(z(:)); 0]);

% get z values to sample the colormaps at
z=zmin+(zmax-zmin).*(0:m-1)'/(m-1);

% setup z values of colormaps
ns=size(sea,1);
zs=zmin-zmin*((0:(ns-1))'/(ns-1));
nl=size(land,1);
zl=zmax*((0:(nl-1))'/(nl-1));

% exclude land or sea if not necessary
if(~zmin)
    % no sea
    z0=zl;
    map=land;
elseif(~zmax)
    % no land
    z0=zs;
    map=sea;
else
    % land & sea
    z0=[zs; zl];
    map=[sea; land];
end

% interpolate new colormap
map=interp1q(z0,map,z);

% fix out-of-bounds due to interpolation
map(map<0)=0;
map(map>1)=1;

% output or adjust colormap?
if(nargout)
    varargout{1}=map;
    varargout{2}=z;
else
    caxis(ax,[zmin zmax]);
    colormap(ax,map);
end

end
