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
%              mantlemap(...,'titletype',value,...)
%              mantlemap(...,'mmap_opt1',mmap_val1,...)
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
%     MANTLEMAP(...,'DVRNG',DVLIMITS,...) sets the coloring limits.
%     Anything outside of this range is set to either the maximum or
%     minimum colors in the color map.  The default is set to the limits of
%     the data.  The units are in % dlnv.
%
%     MANTLEMAP(...,'COLORMAP',CMAP,...) alters the colormap used in the
%     map.  The default colormap is 'seis'.  The colormap may be a Nx3 RGB
%     triplet array or a string that may be evaluated to a Nx3 RGB triplet.
%
%     MANTLEMAP(...,'SHOWCOLORBAR',LOGICAL,...) turns on/off the drawing of
%     a colorbar.  The default is TRUE.
%
%     MANTLEMAP(...,'TITLETYPE',VALUE,...) changes the automatic title
%     info.  Valid values are 0 (or false), 1 (or true -- the default), 2,
%     or 3.  False shows no title.  The others:
%      1: Model Depth   or   Model Min(Depth)-Max(Depth)
%      2: Min(Depth)-Max(Depth)
%      3: Mean(Depths)
%
%     MANTLEMAP(...,'MMAP_OPT1',MMAP_VAL1,...) controls the mapping
%     using MMAP options.  See MMAP or the Examples section below.
%
%     AX=MANTLEMAP(...) returns the axes handle for the map.
%
%    Notes:
%
%    Examples:
%     % Turn off tick labels:
%     mantlemap('dep',1000,'go',{'xticklabel',[],'yticklabel',[]});
%
%     % Compare the CMB region for 4 different P-wave models:
%     figure('color','w');
%     model={'dz04' 'pri05p' 'hmsl06p' 'mitp08'};
%     for i=1:4
%         ax=subplot(2,2,i);
%         mantlemap('dvrng',[-2 2],'mo',model{i},...
%             'parent',ax,'cb',false,'fg','k');
%     end
%
%    See also: MANTLEDV, AVAILABLE_3DMODELS, MANTLEPROFILE

%     Version History:
%        Aug.  4, 2010 - initial version
%        Feb. 24, 2011 - slightly better axes handle code, better labels
%        Aug. 26, 2013 - use mmap image option (big code reduction)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 26, 2013 at 20:30 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:mantlemap:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'m' 'S20RTS' 'd' 2850 'lar' [-90 90] 'lor' [0 360] ...
    'dst' 100 'last' 1 'lost' 1 'dvrng' [] 'cmap' 'seis' ...
    'cb' true 'tt' true} varargin];

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:mantlemap:badOption',...
        'All Options must be specified with a string!');
end

% check option/value pairs
delete=false(1,numel(varargin));
for i=1:2:numel(varargin)
    % skip empty by default (but still checking option exists)
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end
    val=varargin{i+1};
    
    % check option is available
    switch lower(varargin{i})
        case {'model' 'mo' 'm'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(ischar(val))
                model=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'MODEL must be a string!');
            end
        case {'depth' 'dep' 'd'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>=0)
                deprng=[val val];
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTH must be a positive scalar!');
            end
        case {'depthrange' 'deprng' 'dr'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]))
                deprng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTHRANGE must be [depth_lo depth_hi]!');
            end
        case {'latituderange' 'latrng' 'lar'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]) ...
                    && all(abs(val)<=90))
                latrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LATRANGE must be [lat_lo lat_hi]!');
            end
        case {'longituderange' 'lonrng' 'lor'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]))
                lonrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LONRANGE must be [lon_lo lon_hi]!');
            end
        case {'depthstep' 'depstep' 'depinc' 'dst'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                depstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DEPTHSTEP must be a positive scalar!');
            end
        case {'latitudestep' 'latstep' 'latinc' 'last'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                latstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LATSTEP must be a positive scalar!');
            end
        case {'longitudestep' 'lonstep' 'loninc' 'lost'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                lonstep=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'LONSTEP must be a positive scalar!');
            end
        case {'clim' 'dvrng'}
            delete(i:i+1)=true;
            if(skip)
                dvrng=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                dvrng=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'DVRNG must be [%%dv_lo %%dv_hi]!');
            end
        case {'colormap' 'cmap'}
            delete(i:i+1)=true;
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
        case {'showcolorbar' 'colorbar' 'cb'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                showcb=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'SHOWCOLORBAR must be TRUE or FALSE!');
            end
        case {'titletype' 'tt'}
            delete(i:i+1)=true;
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                titletype=val;
            else
                error('seizmo:mantlemap:badInput',...
                    'TITLETYPE must be a scalar value!');
            end
    end
end
varargin(delete)=[];

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

% a couple mmap default changes
% - use min/max of lat/lon as the map boundary
% - do not show land/ocean
varargin=[{'po' {'lat' [min(lat(:)) max(lat(:))] ...
    'lon' [min(lon(:)) max(lon(:))]} 'l' false 'o' false} varargin];

% draw map
ax=mmap('image',{lat lon dv},varargin{:});

% pretty it up
colormap(ax,cmap);
if(~isempty(dvrng)); set(ax,'clim',dvrng); end

% extract foreground color
fg=get(findobj(ax,'tag','m_grid_box'),'color');

% colorbar & title
if(showcb)
    c=colorbar('southoutside','peer',ax,'xcolor',fg,'ycolor',fg);
    [wtype,wtype]=mantledv(model,lat(1),lon(1),dep(1)); % wave type
    xlabel(c,['% \delta{}lnv_' wtype],'color',fg);
end
switch titletype
    case 1
        if(diff(deprng))
            title(ax,[upper(model) '  ' num2str(deprng(1)) '-' ...
                num2str(deprng(2)) 'km'],'color',fg);
        else
            title(ax,[upper(model) '  ' num2str(deprng(1)) 'km'],...
                'color',fg);
        end
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
