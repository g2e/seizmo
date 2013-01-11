function [varargout]=map_cmb_profiles(pf,field,clim,cmap,varargin)
%MAP_CMB_PROFILES    Map CMB profile measurements
%
%    Usage:    map_cmb_profiles(pf,field)
%              map_cmb_profiles(pf,field,clim)
%              map_cmb_profiles(pf,field,clim,cmap)
%              map_cmb_profiles(pf,field,clim,cmap,...,'option',value,...)
%              ax=map_cmb_profiles(...)
%
%    Description:
%     MAP_CMB_PROFILES(PF,FIELD) plots the profiles in struct PF from
%     SLOWDECAYPAIRS or SLOWDECAYPROFILES, colored by the values given in
%     the field FIELD.  FIELD may be 'SLOW', 'CSLOW', 'DECAY' or 'CDECAY'.
%
%     MAP_CMB_PROFILES(PF,FIELD,CLIM) adjusts the color limits for the
%     profile coloring.  The default is [] (empty matrix) and will use the
%     data limits.  CLIM should be specified as [MIN MAX].
%
%     MAP_CMB_PROFILES(PF,FIELD,CLIM,CMAP) specifies the colormap used in
%     coloring the profile segments.  The default is to use the colormap of
%     the figure plotted in.  CMAP is expected to be an Nx3 rgb color array
%     as output by most colormap functions (eg. JET).
%
%     MAP_CMB_PROFILES(PF,FIELD,CLIM,CMAP,...,'OPTION',VALUE,...) allows
%     passing options for mapping.  See MMAP for details.
%
%     AX=MAP_CMB_PROFILES(...) returns the axis handle AX of the plot.
%
%    Notes:
%
%    Examples:
%     % Plot corrected profile slowness with color limits of 4.4-4.7:
%     map_cmb_profiles(pf,'cslow',[4.4 4.7])
%
%    See also: MMAP, SLOWDECAYPAIRS, SLOWDECAYPROFILES,
%              PLOT_CMB_PDF, PLOT_CMB_MEASUREMENTS

%     Version History:
%        Dec. 12, 2010 - initial version
%        Dec. 13, 2010 - fix wrap-around issue
%        Jan. 18, 2011 - new pf fields, fix single profile bug
%        Feb.  1, 2011 - updated See also section, fix coloring of
%                        colorbar, label colorbar, add title
%        Feb. 10, 2011 - update for maplocations=>mmap
%        Mar. 30, 2011 - remove sign flip of decay constant, doc updates
%        Oct. 10, 2012 - bugfix: azimuthal mean via azmean
%        Oct. 11, 2012 - drop corrections field requirement
%        Oct. 16, 2012 - use 27deg over 44deg backprojection to cmb as this
%                        is where paths differ the most (they are near
%                        parallel for the first 20deg along the cmb)
%        Nov. 21, 2012 - no yticks bugfix, colormap input
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 21, 2012 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));
if(nargin>4 && mod(nargin,2))
    error('seizmo:map_cmb_profiles:badNumInputs',...
        'Unpaired Option/Value!');
end

% check profiles struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','corrcoef','freq','phase','runname','dirname',...
    'time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:map_cmb_profiles:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check field
if(~ischar(field) || ...
        ~any(strcmpi(field,{'slow' 'cslow' 'decay' 'cdecay'})))
    error('seizmo:map_cmb_profiles:badInput',...
        'FIELD must be ''SLOW'', ''CSLOW'', ''DECAY'', or ''CDECAY''!');
end
field=lower(field);

% look out for multiple phases
phases=unique({pf.phase}');
if(~isscalar(phases))
    error('seizmo:map_cmb_profiles:badInput',...
        'Multiple phases in dataset not allowed!');
end

% check clim & cmap
if(nargin>2 && ~isempty(clim) && (~isequal(size(clim),[1 2]) ...
        || ~isreal(clim) || ~isnumeric(clim)))
    error('seizmo:map_cmb_profiles:badInput',...
        'CLIM input incorrect! Must be [ZMIN ZMAX]!');
elseif(nargin>3 && ~isempty(cmap) && (size(cmap,2)~=3 ...
        || ~isreal(cmap) || ~isnumeric(cmap) ...
        || any(cmap(:)<0 | cmap(:)>1)))
    error('seizmo:map_cmb_profiles:badInput',...
        'CMAP input incorrect! Must be a Nx3 array of values from 0-1!');
end

% unique stations & events
st=unique(cell2mat({pf.st}'),'rows');
ev=unique(cell2mat({pf.ev}'),'rows');

% call mmap
ax=mmap('stations',st(:,1:2),'events',ev(:,1:2),varargin{:});

% slowness to color
slow=[pf.(field)].';
if(nargin<3 || isempty(clim))
    if(numel(pf)>1)
        clim=[min(slow) max(slow)];
    else
        % only fails if zero
        clim=[slow*.9 slow*1.1];
    end
end
if(nargin<4 || isempty(cmap)); cmap=colormap(ax); end
color=z2c(slow,cmap,clim);

% loop over profiles, plotting each
hold(ax,'on');
for a=1:numel(pf)
    % get profile endpoints
    evla=pf(a).ev(1,1);
    evlo=pf(a).ev(1,2);
    distmin=min(pf(a).delaz(:,4))-27*111.19;
    distmax=max(pf(a).delaz(:,4))-27*111.19;
    azavg=azmean(pf(a).delaz(:,2));
    [lat1,lon1]=vincentyfwd(evla,evlo,distmin,azavg);
    [lat2,lon2]=vincentyfwd(evla,evlo,distmax,azavg);
    
    % plot profile
    [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2,10);
    lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
    m_line(lon',lat','linewi',2,'color',color(a,:));
end
hold(ax,'off');

% add title
tc=get(findobj(ax,'tag','m_grid_box'),'color');
if(iscell(tc)); tc=tc{1}; end
switch lower(field)
    case 'slow'
        title(ax,[phases{1} ' Ray Parameter'],'color',tc);
    case 'cslow'
        title(ax,[phases{1} ' Corrected Ray Parameter'],'color',tc);
    case 'decay'
        title(ax,[phases{1} ' Decay Constant'],'color',tc);
    case 'cdecay'
        title(ax,[phases{1} ' Corrected Decay Constant'],'color',tc);
end

% add colormap
colormap(ax,cmap);
set(ax,'clim',clim);
cb=colorbar('peer',ax);
set(cb,'xcolor',tc,'ycolor',tc);
switch lower(field)
    case {'slow' 'cslow'}
        set(get(cb,'xlabel'),'string','s/^o');
end

% move stations/events forward
sta=findobj(ax,'tag','stations');
evt=findobj(ax,'tag','events');
movekids(evt,'front');
movekids(sta,'front');

% optional output
if(nargout); varargout{1}=ax; end

end
