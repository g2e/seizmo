function [varargout]=mantlemap(varargin)
%MANTLEMAP    3D mantle model map (aka depth slice)
%
%    Usage:    mantlemap()
%              mantlemap(...,'model',model,...)
%              mantlemap(...,'depth',depth,...)
%              mantlemap(...,'deprng',deplimits,...)
%              mantlemap(...,'latrng',latlimits,...)
%              mantlemap(...,'lonrng',lonlimits,...)
%              mantlemap(...,'depinc',increment,...)
%              mantlemap(...,'latinc',increment,...)
%              mantlemap(...,'loninc',increment,...)
%              mantlemap(...,'dvrng',dvlimits,...)
%              mantlemap(...,'colormap',cmap,...)
%              mantlemap(...,'gshhs',res,...)
%              mantlemap(...,'proj',proj,...)
%              mantlemap(...,'projopt',{'opt',val,...},...)
%              mantlemap(...,'gridopt',{'opt',val,...},...)
%              mantlemap(...,'fgcolor',color,...)
%              mantlemap(...,'bgcolor',color,...)
%              mantlemap(...,'coast',color,...)
%              mantlemap(...,'border',color,...)
%              mantlemap(...,'showcolorbar',logical,...)
%              mantlemap(...,'axis',ax,...)
%              ax=mantlemap(...)
%
%    Description:
%     MANTLEMAP() plots S20RTS at 2850km depth for the entire globe in a
%     Robinson projection map.  The model is sampled at 1deg increments and
%     the colormap is set to SEIS.  Coastlines are shown in black.  The map
%     is drawn in a new figure with a white axis on a black background.  A
%     colorbar is drawn below the map and is in units of % dlnv.
%
%     MANTLEMAP(...,'MODEL',MODEL,...) indicates the 3D mantle model for
%     this map.  The default is 'S20RTS'.  Call AVAILABLE_3DMODELS for a
%     list.
%
%     MANTLEMAP(...,'DEPTH',DEPTH,...) specifies a specific depth to sample
%     the mantle model at.  This will override the more general DEPRNG
%     option below.  Units are in kilometers.
%
%     MANTLEMAP(...,'DEPRNG',DEPLIMITS,...) specifies the range of the
%     depth values to sample the mantle model at (the depths in this range
%     will be averaged).  The default is [2850 2850].  So the mantle model
%     is only sampled at a single depth by default.  If a range is
%     specified the model is sampled every 100km unless altered by the
%     DEPINC option.  Units are in kilometers.
%
%     MANTLEMAP(...,'LATRNG',LATLIMITS,...) specifies the range of the
%     latitude values to sample the mantle model at.  The default is
%     [-90 90].  The model is sampled every 1 degree unless altered by the
%     LATINC option.
%
%     MANTLEMAP(...,'LONRNG',LONLIMITS,...) specifies the range of the
%     longitude values to sample the mantle model at.  The default is
%     [0 360].  The model is sampled every 1 degree unless altered by the
%     LONINC option.
%
%     MANTLEMAP(...,'DEPINC',INCREMENT,...) sets the depth increment
%     for sampling the mantle model for creating an average.  The default
%     is 100 km.
%
%     MANTLEMAP(...,'LATINC',INCREMENT,...) sets the latitude increment
%     for sampling the mantle model.  The default is 1 degree.
%
%     MANTLEMAP(...,'LONINC',INCREMENT,...) sets the longitude increment
%     for sampling the mantle model.  The default is 1 degree.
%
%     MANTLEMAP(...,'COLORMAP',CMAP,...) alters the colormap used in the
%     map.  The default colormap is 'seis'.  The colormap may be a Nx3 RGB
%     triplet array or a string that may be evaluated to a Nx3 RGB triplet.
%
%     MANTLEMAP(...,'CLIM',DVLIMITS,...) sets the coloring limits.
%     Anything outside of this range is set to either the maximum or
%     minimum colors in the color map.  The default is set to the limits of
%     the data.  The units are in % dlnv.
%
%     MANTLEMAP(...,'GSHHS',RES,...) sets the GSHHS coastline and
%     political boundaries resolution.  The values can be 'c', 'l', 'i',
%     'h', 'f', or 'o' (for 'off').  The default GSHHS resolution is 'o'
%     (off) which calls M_COAST and does not draw political borders.
%
%     MANTLEMAP(...,'PROJ',PROJ,...) defines the map projection.  See
%     M_PROJ('SET') for possible projections.  The default PROJ is
%     'Robinson'.
%
%     MANTLEMAP(...,'PROJOPT',{'OPT',VAL,...},...) passes additional
%     options to M_PROJ (like the lat/lon boundaries of the map).  The
%     options must be wrapped in a cell array!  See M_PROJ('get',PROJ) for
%     a list of possible options for the set projection (see 'PROJ' option
%     for the default projection and altering it).  The default will create
%     a map to the limits of the lat/lon ranges of the dv data.
%
%     MANTLEMAP(...,'GRIDOPT',{'OPT',VAL,...},...) passes options to
%     M_GRID (like the lat/lon ticks of the map, etc).  The options must be
%     wrapped in a cell array!  See M_GRID('get') for a list of possible
%     options and M_GRID('set') for their defaults.  The default is no
%     options.
%
%     MANTLEMAP(...,'FGCOLOR',COLOR,...) specifies the foreground color
%     of the map.  The default is 'w'.  If BGCOLOR is specified and FGCOLOR
%     is not, then FGCOLOR will be set using INVERTCOLOR.
%
%     MANTLEMAP(...,'BGCOLOR',COLOR,...) specifies the background color
%     of the map.  The default is 'k'.  If FGCOLOR is specified and BGCOLOR
%     is not, then BGCOLOR will be set using INVERTCOLOR.
%
%     MANTLEMAP(...,'COAST',COLOR,...) specifies the color of the coastline
%     in the map.  The default is 'k' (black).
%
%     MANTLEMAP(...,'BORDER',COLOR,...) specifies the color of the
%     political borders in the map.  The default is [.5 0 0].
%
%     MANTLEMAP(...,'SHOWCOLORBAR',LOGICAL,...) turns on/off the drawing of
%     a colorbar.  The default is TRUE.
%
%     MANTLEMAP(...,'AXIS',AX,...) sets the axes to draw in.  This is
%     useful for subplots, guis, etc.  The default draws the map in a new
%     figure.
%
%     AX=MANTLEMAP(...) returns the axes handle for the map.
%
%    Notes:
%
%    Examples:
%     % compare the CMB region for 4 different P-wave models
%      figure('color','w');
%      model={'dz04' 'pri05p' 'hmsl06p' 'mitp08'};
%      for i=1:4
%        ax=subplot(2,2,i);
%        mantlemap('clim',[-2 2],'mo',model{i},...
%          'ax',ax,'cb',false,'fg','w');
%      end
%
%    See also: MANTLEDV, AVAILABLE_3DMODELS

%     Version History:
%        Aug.  4, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  4, 2010 at 20:30 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:mantlemap:badNumInputs',...
        'Unpaired Option/Value!');
end

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% option defaults
varargin=[{'m' 'S20RTS' 'd' 2850 'lar' [-90 90] 'lor' [0 360] ...
    'dst' 100 'last' 1 'lost' 1 'dvrng' [] 'cmap' 'seis' 'g' 'o' ...
    'proj' 'robinson' 'po' [] 'go' [] 'fg' [] 'bg' [] 'c' 'k' ...
    'b' [.5 0 0] 'cb' true 'tt' true 'a' []} varargin];

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:mantlemap:badOption',...
        'All Options must be specified with a string!');
end

% check option/value pairs
for i=1:2:numel(varargin)
    % skip empty by default (but still checking option exists)
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end
    val=varargin{i+1};
    
    % check option is available
    switch lower(varargin{i})
        case {'model' 'mo' 'm'}
            if(skip); continue; end
            if(ischar(val))
                model=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'MODEL must be a string!');
            end
        case {'depth' 'dep' 'd'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>=0)
                deprng=[val val];
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTH must be a positive scalar!');
            end
        case {'depthrange' 'deprng' 'dr'}
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]))
                deprng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTHRANGE must be [depth_lo depth_hi]!');
            end
        case {'latituderange' 'latrng' 'lar'}
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]) ...
                    && all(abs(val)<=90))
                latrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LATRANGE must be [lat_lo lat_hi]!');
            end
        case {'longituderange' 'lonrng' 'lor'}
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]))
                lonrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LONRANGE must be [lon_lo lon_hi]!');
            end
        case {'depthstep' 'depstep' 'depinc' 'dst'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                depstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTHSTEP must be a positive scalar!');
            end
        case {'latitudestep' 'latstep' 'latinc' 'last'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                latstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LATSTEP must be a positive scalar!');
            end
        case {'longitudestep' 'lonstep' 'loninc' 'lost'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                lonstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LONSTEP must be a positive scalar!');
            end
        case {'clim' 'dvrng'}
            if(skip)
                dvrng=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                dvrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DVRANGE must be [%%dv_lo %%dv_hi]!');
            end
        case {'colormap' 'cmap'}
            if(skip); continue; end
            if(isreal(val) && ndims(val)==2 ...
                    && size(val,2)==3 ...
                    && all(val(:)>=0 & val(:)<=1))
                cmap=val;
            elseif(ischar(val) && isvector(val) && size(val,1)==1)
                cmap=val;
            else
                error('seizmo:mantlemap:badInput',...
                    ['COLORMAP must be a colormap function\n'...
                    'string or a Nx3 RGB triplet array!']);
            end
        case {'gshhs' 'g' 'resolution' 'res'}
            if(skip); continue; end
            if(ischar(val) && numel(val)==1 ...
                    && any(strcmpi(val,{'o' 'c' 'l' 'i' 'h' 'f'})))
                gshhs=lower(val);
            else
                error('seizmo:mantlemap:badInput',...
                    'GSHHS option must be c, l, i, h, or f !');
            end
        case {'projection' 'proj' 'p'}
            if(skip); continue; end
            if(ischar(val) && ndims(val)==2 && size(val,1)==1)
                proj=lower(val);
            else
                error('seizmo:mantlemap:badInput',...
                    'PROJECTION option must be a string!');
            end
        case {'projopt' 'popt' 'po'}
            if(skip)
                popt={};
            elseif(iscell(val) && iscellstr(val(1:2:end)))
                popt=val;
            else
                error('seizmo:mantlemap:badInput',...
                    ['PROJOPT option must be a cell array of ' ...
                    '''option''/value pairs!']);
            end
        case {'gridopt' 'gopt' 'go'}
            if(skip)
                gopt={};
            elseif(iscell(val) && iscellstr(val(1:2:end)))
                gopt=val;
            else
                error('seizmo:mantlemap:badInput',...
                    ['GRIDOPT option must be a cell array of ' ...
                    '''option''/value pairs!']);
            end
        case {'fgcolor' 'fg'}
            if(skip)
                fg=[];
            elseif(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                fg=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'FGCOLOR must be a colorname or RGB triplet!');
            end
        case {'bgcolor' 'bg'}
            if(skip)
                bg=[];
            elseif(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                bg=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'BGCOLOR must be a colorname or RGB triplet!');
            end
        case {'coastcolor' 'coast' 'c'}
            if(skip)
                coast=[];
            elseif(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                coast=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'COASTCOLOR must be a colorname or RGB triplet!');
            end
        case {'bordercolor' 'border' 'b'}
            if(skip); continue; end
            if(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                border=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'BORDERCOLOR must be a colorname or RGB triplet!');
            end
        case {'showcolorbar' 'colorbar' 'cb'}
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                showcb=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'SHOWCOLORBAR must be TRUE or FALSE!');
            end
        case {'titletype' 'tt'}
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                titletype=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'TITLETYPE must be a scalar value!');
            end
        case {'axis' 'ax' 'a'}
            if(skip)
                ax=[];
            else
                ax=val;
            end
        otherwise
            error('seizmo:mantlemap:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% fix fg/bg colors
if(isempty(fg))
    if(isempty(bg))
        fg='w'; bg='k';
    else
        fg=invertcolor(bg,true);
    end
elseif(isempty(bg))
    bg=invertcolor(fg,true);
end
if(isempty(coast))
    coast=fg;
end

% convert colornames
if(ischar(fg)); fg=name2rgb(fg); end
if(ischar(bg)); bg=name2rgb(bg); end
if(ischar(coast)); coast=name2rgb(coast); end
if(ischar(border)); border=name2rgb(border); end

% lat/lon grid
lat0=latrng(1):latstep:latrng(2);
lon0=lonrng(1):lonstep:lonrng(2);
[lon,lat]=meshgrid(lon0,lat0);

% loop over depths, get average %dv
dep=deprng(1):depstep:deprng(2);
dv=zeros(size(lat));
for i=1:numel(dep)
    dv=dv+mantledv(model,lat,lon,dep(i));
end
dv=dv./numel(dep).*100;

% image to pcolor
lat=lat-latstep/2;
lon=lon-lonstep/2;
lat=lat([1:end end],[1:end end]);
lon=lon([1:end end],[1:end end]);
lat(end,:)=lat(end,:)+latstep;
lon(:,end)=lon(:,end)+lonstep;
dv=dv([1:end end],[1:end end]);

% force >90 lat to 90
lat(lat<-90)=-90;
lat(lat>90)=90;

% default projection options
if(isempty(popt))
    % use min/max of lat/lon
    minlat=min(lat(:));
    maxlat=max(lat(:));
    minlon=min(lon(:));
    maxlon=max(lon(:));
    popt={'lat',[minlat maxlat],'lon',[minlon maxlon]};
end

% setup axis
if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    % new figure
    figure('color',bg);
    ax=gca;
else
    axes(ax);
    h=get(ax,'children'); delete(h);
    h=findobj(get(get(ax,'parent'),'children'),'peer',ax); delete(h);
end

% setup projection
m_proj(proj,popt{:});
set(ax,'color',bg);

% plot mantle map
hold on
if(any(lon(:)>MAP_VAR_LIST.longs(1) ...
        & lon(:)<MAP_VAR_LIST.longs(2)))
    m_pcolor(lon,lat,dv);
end
if(any(lon(:)-360>MAP_VAR_LIST.longs(1) ...
        & lon(:)-360<MAP_VAR_LIST.longs(2)))
    m_pcolor(lon-360,lat,dv);
end
if(any(lon(:)+360>MAP_VAR_LIST.longs(1) ...
        & lon(:)+360<MAP_VAR_LIST.longs(2)))
    m_pcolor(lon+360,lat,dv);
end

% pretty it up
shading flat;
colormap(cmap);
if(~isempty(dvrng)); set(ax,'clim',dvrng); end
hold off

% plot coasts & borders
if(strcmpi(gshhs,'o'))
    m_coast('color',coast);
else
    m_gshhs([gshhs 'c'],'color',coast);
    m_gshhs([gshhs 'b'],'color',border);
end
m_grid('color',fg,gopt{:});

% hackery to color the map if the above fails
set(findobj(ax,'tag','m_grid_color'),'facecolor',bg);

% colorbar & title
if(showcb)
    c=colorbar('southoutside','peer',ax,'xcolor',fg,'ycolor',fg);
    xlabel(c,'% dlnv','color',fg);
end
switch titletype
    case 1
        title(ax,[upper(model) '  ' num2str(deprng(1)) '-' ...
            num2str(deprng(2)) 'km'],'color',fg);
    case 2
        title(ax,[num2str(deprng(1)) '-' num2str(deprng(2)) 'km'],...
            'color',fg);
    case 3
        title(ax,[num2str(mean(dep)) 'km'],'color',fg);
end

% export basic info
userdata.model=model;
userdata.depths=dep;
set(ax,'userdata',userdata);

% return figure handle
set(ax,'tag','mantlemap');
if(nargout); varargout{1}=ax; end

end
