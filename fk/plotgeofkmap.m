function [varargout]=plotgeofkmap(map,popt,dblim,zerodb,fgcolor,bgcolor,ax)
%PLOTGEOFKMAP    Plots frequency-slowness-position response info
%
%    Usage:    plotgeofkmap(map)
%              plotgeofkmap(map,projopt)
%              plotgeofkmap(map,projopt,dblim)
%              plotgeofkmap(map,projopt,dblim,zerodb)
%              plotgeofkmap(map,projopt,dblim,zerodb,fgcolor,bgcolor)
%              plotgeofkmap(map,projopt,dblim,zerodb,fgcolor,bgcolor,ax)
%              ax=plotgeofkmap(...)
%
%    Description: PLOTGEOFKMAP(MAP) plots the frequency-slowness-position
%     response in geofk struct MAP.  See a geofk function like
%     GEOFKXCVOLUME for details on the struct.  The response is plotted on
%     a map with a Hammer-Aitoff projection and the map limits are scaled
%     to fit the response & station positions.  Note that the response
%     positions should form a regular grid (using a function like
%     MESHGRID).  This plots GSHHS coastlines and borders in low-resolution
%     which may take a few moments - please be patient.
%
%     PLOTGEOFKMAP(MAP,PROJOPT) allows passing options to M_PROJ.  See
%     M_PROJ('SET') for possible projections and See M_PROJ('GET',PROJ) for
%     a list of possible additional options specific to that projection.
%
%     PLOTGEOFKMAP(MAP,PROJOPT,DBLIM) sets the dB limits for coloring the
%     response info.  The default is [-12 0] for the default ZERODB (see
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTGEOFKMAP(MAP,PROJOPT,DBLIM,ZERODB) changes what 0dB corresponds
%     to in the plot.  The allowed values are 'min', 'max', 'median', &
%     'abs'.  The default is 'max'.
%
%     PLOTGEOFKMAP(MAP,PROJOPT,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies
%     foreground and background colors of the plot.  The default is 'w' for
%     FGCOLOR & 'k' for BGCOLOR.  Note that if one is specified and the
%     other is not, an opposing color is found using INVERTCOLOR.  The
%     color scale is also changed so the noise clip is at BGCOLOR.
%
%     PLOTGEOFKMAP(MAP,PROJOPT,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the
%     axes to draw in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTGEOFKMAP(...)
%
%    Notes:
%
%    Examples:
%     In search of the 26s microseism:
%      [lat,lon]=meshgrid(-60:60,-60:60);
%      zgeo=geofkxcvolume(xcdata,[lat(:) lon(:)],27:33,[1/26.3 1/26]);
%      zgeo0=geofkvol2map(zgeo);
%      plotgeofkmap(zgeo0);
%
%    See also: GEOFKFREQSLIDE, GEOFKSLOWSLIDE, GEOFKVOL2MAP, GEOFKXCVOLUME,
%              GEOFKXCHORZVOLUME, CHKGEOFKSTRUCT

%     Version History:
%        June 25, 2010 - initial version
%        July  1, 2010 - no land cover
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July  1, 2010 at 15:45 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check fk struct
error(chkgeofkstruct(map));

% don't allow array/volume
if(~isscalar(map) || any(map.volume))
    error('seizmo:plotgeofkmap:badInput',...
        'MAP must be a scalar geofk struct and not a volume!');
end

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotgeofkmap:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<3 || isempty(dblim))
    switch zerodb
        case {'min' 'median'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotgeofkmap:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale response
switch zerodb
    case 'min'
        map.normdb=map.normdb+min(map.response(:));
        map.response=map.response-min(map.response(:));
    case 'max'
        map.normdb=map.normdb+max(map.response(:));
        map.response=map.response-max(map.response(:));
    case 'median'
        map.normdb=map.normdb+median(map.response(:));
        map.response=map.response-median(map.response(:));
    case 'abs'
        map.response=map.response+map.normdb;
        map.normdb=0;
end

% check colors
if(nargin<5); fgcolor='w'; bgcolor='k'; end
if(nargin<6)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end
if(isempty(fgcolor))
    if(isempty(bgcolor))
        fgcolor='w'; bgcolor='k';
    else
        fgcolor=invertcolor(bgcolor,true);
    end
elseif(isempty(bgcolor))
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(nargin<7 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
    h=get(ax,'children'); delete(h);
end

% map colors & coast/border res
gshhs='c';
%ocean=[0.3 0.6 1];
%land=[0.4 0.6 0.2];
%border=[0.5 0 0];
ocean=bgcolor;
land=fgcolor;
border=fgcolor;

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% reshape response & account for pcolor
nlat=numel(unique(map.latlon(:,1)));
nlon=numel(unique(map.latlon(:,2)));
map.latlon=reshape(map.latlon,[nlon nlat 2]);
map.response=reshape(map.response,[nlon nlat]);
latstep=map.latlon(1,2,1)-map.latlon(1,1,1);
lonstep=map.latlon(2,1,2)-map.latlon(1,1,2);
map.latlon(:,:,1)=map.latlon(:,:,1)-latstep/2;
map.latlon(:,:,2)=map.latlon(:,:,2)-lonstep/2;
map.latlon=map.latlon([1:end end],[1:end end],:);
map.latlon(:,end,1)=map.latlon(:,end,1)+latstep;
map.latlon(end,:,2)=map.latlon(end,:,2)+lonstep;
map.response=map.response([1:end end],[1:end end]);

% get max/min lat/lon of map & stations
minlat=min([map.stla; min(map.latlon(:,:,1))']);
maxlat=max([map.stla; max(map.latlon(:,:,1))']);
minlon=min([map.stlo; min(map.latlon(:,:,2))']);
maxlon=max([map.stlo; max(map.latlon(:,:,2))']);

% default/check projopt
if(nargin<2 || isempty(popt))
    popt={'hammer','lat',[minlat maxlat],'lon',[minlon maxlon]};
end
if(ischar(popt)); popt=cellstr(popt); end
if(~iscell(popt))
    error('seizmo:plotgeofkmap:badInput',...
        'PROJOPT must be a cell array of args for M_PROJ!');
end

% setup projection
m_proj(popt{:});
set(ax,'color',ocean);

% plot geofk response
hold on
if(any(map.latlon(:,:,2)>MAP_VAR_LIST.longs(1) ...
        & map.latlon(:,:,2)<MAP_VAR_LIST.longs(2)))
    m_pcolor(map.latlon(:,:,2),map.latlon(:,:,1),double(map.response));
end
if(any(map.latlon(:,:,2)-360>MAP_VAR_LIST.longs(1) ...
        & map.latlon(:,:,2)-360<MAP_VAR_LIST.longs(2)))
    m_pcolor(map.latlon(:,:,2)-360,map.latlon(:,:,1),double(map.response));
end
if(any(map.latlon(:,:,2)+360>MAP_VAR_LIST.longs(1) ...
        & map.latlon(:,:,2)+360<MAP_VAR_LIST.longs(2)))
    m_pcolor(map.latlon(:,:,2)+360,map.latlon(:,:,1),double(map.response));
end

% modify
shading flat;
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(hsvcustom(hsv));
end
set(ax,'clim',dblim);
hold off

% now add coastlines and political boundaries
%m_gshhs([gshhs 'c'],'patch',land);
%m_gshhs([gshhs 'b'],'color',border);
m_gshhs([gshhs 'c'],'color',land);
m_gshhs([gshhs 'b'],'color',border);
m_grid('color',fgcolor);

% hackery to color oceans at large when the above fails
set(findobj(ax,'tag','m_grid_color'),'facecolor',ocean);

% wrap station longitudes to within 180deg of plot center
while(any(abs(map.stlo-mean(MAP_VAR_LIST.longs))>180))
    map.stlo(map.stlo<MAP_VAR_LIST.longs(1))=...
        map.stlo(map.stlo<MAP_VAR_LIST.longs(1))+360;
    map.stlo(map.stlo>MAP_VAR_LIST.longs(2))=...
        map.stlo(map.stlo>MAP_VAR_LIST.longs(2))-360;
end

% add stations
hold on
h=m_scatter(map.stlo,map.stla,[],'y','filled',...
    'markeredgecolor','k');
set(h,'tag','stations');
hold off

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fgcolor,'ycolor',fgcolor);
xlabel(c,'dB','color',fgcolor);
fmin=min(map.freq); fmax=max(map.freq);
smn=min(map.horzslow); smx=max(map.horzslow);
title(ax,{[] ['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period    :  ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Horiz. Slowness :  ' num2str(smn) ' to ' num2str(smx) ' s/\circ'] ...
    ['0 dB = ' num2str(map.normdb) 'dB'] []},'color',fgcolor);

% set zerodb & dblim in userdata
% - this is for updategeofkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata);

% return figure handle
if(nargout); varargout{1}=ax; end

end
