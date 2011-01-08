function [varargout]=plotgeofkarf(arf,popt,dblim,zerodb,fgcolor,bgcolor,ax)
%PLOTGEOFKARF    Plots geofk array response
%
%    Usage:    plotgeofkarf(arf)
%              plotgeofkarf(arf,projopt)
%              plotgeofkarf(arf,projopt,dblim)
%              plotgeofkarf(arf,projopt,dblim,zerodb)
%              plotgeofkarf(arf,projopt,dblim,zerodb,fgcolor,bgcolor)
%              plotgeofkarf(arf,projopt,dblim,zerodb,fgcolor,bgcolor,ax)
%              ax=plotgeofkarf(...)
%
%    Description: PLOTGEOFKARF(ARF) plots the frequency-slowness-position
%     beam data in geofk struct ARF.  See a geofk function like
%     GEOFKXCVOLUME for details on the struct.  The data is plotted on a
%     map with a Hammer-Aitoff projection and the map limits are scaled to
%     fit the beam data & station positions.  Note that the beam data
%     positions should form a regular grid (using a function like
%     MESHGRID).  This plots GSHHS coastlines and borders in low-resolution
%     which may take a few moments - please be patient.
%
%     PLOTGEOFKARF(ARF,PROJOPT) allows passing options to M_PROJ.  See
%     M_PROJ('SET') for possible projections and See M_PROJ('GET',PROJ) for
%     a list of possible additional options specific to that projection.
%
%     PLOTGEOFKARF(ARF,PROJOPT,DBLIM) sets the dB limits for coloring the
%     response info.  The default is [-12 0] for the default ZERODB (see
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTGEOFKARF(ARF,PROJOPT,DBLIM,ZERODB) changes what 0dB corresponds
%     to in the plot.  The allowed values are 'min', 'max', 'median', &
%     'abs'.  The default is 'max'.
%
%     PLOTGEOFKARF(ARF,PROJOPT,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies
%     foreground and background colors of the plot.  The default is 'w' for
%     FGCOLOR & 'k' for BGCOLOR.  Note that if one is specified and the
%     other is not, an opposing color is found using INVERTCOLOR.  The
%     color scale is also changed so the noise clip is at BGCOLOR.
%
%     PLOTGEOFKARF(ARF,PROJOPT,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the
%     axes to draw in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTGEOFKARF(...)
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKARF, GEOFKARF2MAP, GEOFKSUBARF, UPDATEGEOFKARF,
%              GEOFKARFSLOWSLIDE, CHKGEOFKARFSTRUCT

%     Version History:
%        July  7, 2010 - update for new struct
%        Oct.  6, 2010 - truncate title if too many ARF locations
%        Oct. 10, 2010 - all plotting functions use proper ax calls, tagged
%                        plots as 'fkmap'
%        Dec.  8, 2010 - use '^o' for deg symbol rather than \circ
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  8, 2010 at 18:25 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf) || any(arf.volume))
    error('seizmo:plotgeofkarf:badInput',...
        'ARF must be a scalar geofk struct and not a volume!');
end

% default/check scaling type
if(nargin<4 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotgeofkarf:badSTYPE',...
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
    error('seizmo:plotgeofkarf:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale response
switch zerodb
    case 'min'
        arf.normdb=arf.normdb+min(arf.beam(:));
        arf.beam=arf.beam-min(arf.beam(:));
    case 'max'
        arf.normdb=arf.normdb+max(arf.beam(:));
        arf.beam=arf.beam-max(arf.beam(:));
    case 'median'
        arf.normdb=arf.normdb+median(arf.beam(:));
        arf.beam=arf.beam-median(arf.beam(:));
    case 'abs'
        arf.beam=arf.beam+arf.normdb;
        arf.normdb=0;
end

% check colors
if(nargin<5);
    fgcolor='w'; bgcolor='k';
elseif(nargin<6)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
else
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
gshhs='l';
%ocean=[0.3 0.6 1];
%land=[0.4 0.6 0.2];
%border=[0.5 0 0];
ocean=bgcolor;
land=fgcolor;
border=fgcolor;

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% reshape beam & account for pcolor
nlat=numel(unique(arf.latlon(:,1)));
nlon=numel(unique(arf.latlon(:,2)));
arf.latlon=reshape(arf.latlon,[nlon nlat 2]);
arf.beam=reshape(arf.beam,[nlon nlat]);
latstep=arf.latlon(1,2,1)-arf.latlon(1,1,1);
lonstep=arf.latlon(2,1,2)-arf.latlon(1,1,2);
arf.latlon(:,:,1)=arf.latlon(:,:,1)-latstep/2;
arf.latlon(:,:,2)=arf.latlon(:,:,2)-lonstep/2;
arf.latlon=arf.latlon([1:end end],[1:end end],:);
arf.latlon(:,end,1)=arf.latlon(:,end,1)+latstep;
arf.latlon(end,:,2)=arf.latlon(end,:,2)+lonstep;
arf.beam=arf.beam([1:end end],[1:end end]);

% get max/min lat/lon of arf & stations
minlat=min([arf.stla; min(arf.latlon(:,:,1))']);
maxlat=max([arf.stla; max(arf.latlon(:,:,1))']);
minlon=min([arf.stlo; min(arf.latlon(:,:,2))']);
maxlon=max([arf.stlo; max(arf.latlon(:,:,2))']);

% default/check projopt
if(nargin<2 || isempty(popt))
    popt={'hammer','lat',[minlat maxlat],'lon',[minlon maxlon]};
end
if(ischar(popt)); popt=cellstr(popt); end
if(~iscell(popt))
    error('seizmo:plotgeofkarf:badInput',...
        'PROJOPT must be a cell array of args for M_PROJ!');
end

% setup projection
axes(ax);
m_proj(popt{:});
set(ax,'color',ocean);

% plot geofk beam
hold(ax,'on');
if(any(arf.latlon(:,:,2)>MAP_VAR_LIST.longs(1) ...
        & arf.latlon(:,:,2)<MAP_VAR_LIST.longs(2)))
    m_pcolor(arf.latlon(:,:,2),arf.latlon(:,:,1),double(arf.beam),...
        'parent',ax);
end
if(any(arf.latlon(:,:,2)-360>MAP_VAR_LIST.longs(1) ...
        & arf.latlon(:,:,2)-360<MAP_VAR_LIST.longs(2)))
    m_pcolor(arf.latlon(:,:,2)-360,arf.latlon(:,:,1),double(arf.beam),...
        'parent',ax);
end
if(any(arf.latlon(:,:,2)+360>MAP_VAR_LIST.longs(1) ...
        & arf.latlon(:,:,2)+360<MAP_VAR_LIST.longs(2)))
    m_pcolor(arf.latlon(:,:,2)+360,arf.latlon(:,:,1),double(arf.beam),...
        'parent',ax);
end

% modify
shading(ax,'flat');
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(ax,fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(ax,hsvcustom(hsv));
end
set(ax,'clim',dblim);
hold(ax,'off');

% now add coastlines and political boundaries
axes(ax);
%m_gshhs([gshhs 'c'],'patch',land);
%m_gshhs([gshhs 'b'],'color',border);
m_gshhs([gshhs 'c'],'color',land);
m_gshhs([gshhs 'b'],'color',border);
m_grid('color',fgcolor);

% hackery to color oceans at large when the above fails
set(findobj(ax,'tag','m_grid_color'),'facecolor',ocean);

% wrap station longitudes to within 180deg of plot center
while(any(abs(arf.stlo-mean(MAP_VAR_LIST.longs))>180))
    arf.stlo(arf.stlo<MAP_VAR_LIST.longs(1))=...
        arf.stlo(arf.stlo<MAP_VAR_LIST.longs(1))+360;
    arf.stlo(arf.stlo>MAP_VAR_LIST.longs(2))=...
        arf.stlo(arf.stlo>MAP_VAR_LIST.longs(2))-360;
end

% add stations
hold(ax,'on');
h=m_scatter(ax,arf.stlo,arf.stla,[],'y','filled',...
    'markeredgecolor','k');
set(h,'tag','stations');
hold(ax,'off');

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fgcolor,'ycolor',fgcolor);
xlabel(c,'dB','color',fgcolor);
smn=min(arf.horzslow); smx=max(arf.horzslow);
if(arf.nsw<=5)
    titstr=cell(arf.nsw,1);
    for i=1:arf.nsw
        titstr{i}=sprintf(['SLOWNESS: %gs/^o, LAT: %g^o, ' ...
            'LON: %g^o, PERIOD: %gs'],arf.horzslow0(i),...
            arf.latlon0(i,1),arf.latlon0(i,2),1/arf.freq0(i));
    end
else
    titstr{1}=[num2str(arf.nsw) ' Locations'];
end
title(ax,[{[]}; 'Array Response Function @ '; titstr; ...
    ['Number of Stations: ' num2str(arf.nsta)]; ...
    ['Horiz. Slowness : ' num2str(smn) ' to ' num2str(smx) ' s/^o']; ...
    ['0 dB = ' num2str(arf.normdb) 'dB']; {[]}],'color',fgcolor);

% set zerodb & dblim in userdata
% - this is for updategeofkarf
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','geofkarf');

% return figure handle
if(nargout); varargout{1}=ax; end

end
