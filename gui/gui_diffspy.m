function []=gui_diffspy(a)

% check table
if(nargin)
    check_table(a);
else
    % make blank table
    a.align=[];
    a.event.centroidlat=0;
    a.event.centroidlon=0;
end

% spawn figure
f=spawn_diffspy(a);

% Assign the GUI a name to appear in the window title.
set(f,'Name','Core Diffracted Profiles Spy Tool')

% Move the GUI to the center of the screen.
movegui(f,'center')

% Make the GUI visible.
set(f,'Visible','on');

end



function check_table(a)
% check if input table of data looks ok

% required main fields
reqmfields={'align' 'event'};
% required align fields
reqafields={'network' 'station' 'stream' 'channel' 'latitude' 'longitude'...
    'distance' 'azimuth' 'snr' 'cluster' 'tt' 'tterr' 'gtt' 'gtterr'...
    'amp' 'amperr' 'rlow' 'rmean' 'rhigh' 'pol'};
% required event fields (from ndk struct, but just need location)
reqefields={'centroidlat' 'centroidlon'};

if(~isstruct(a))
    error('ALIGN not a struct!');
end

fields=fieldnames(a);
if(~isempty(setdiff(reqmfields,fields)))
    error('ALIGN struct does not contain all required fields!');
end
fields=fieldnames(a.align);
if(~isempty(setdiff(reqafields,fields)))
    error('ALIGN struct does not contain all required fields!');
end
fields=fieldnames(a.event);
if(~isempty(setdiff(reqefields,fields)))
    error('ALIGN struct does not contain all required fields!');
end

end



function [f]=spawn_diffspy(a)

% defaults for data and profile selection
option.cluster='';
option.snrmin=1.8;
option.network='';
option.station='';
option.stream='';
option.channel='';
option.region.azmin=290;
option.region.azmax=310;
option.region.dmin=100;
option.region.dmax=140;
option.profile.ddmin=5;
option.profile.ddmax=inf;
option.profile.dazmin=0;
option.profile.dazmax=15;
option.profile.aspectmin=3;
option.profile.aspectmax=inf;

% create figure and then hide it as it is being constructed
f=figure('visible','off','position',[100 100 600 600]);

% normalized units and coloring
set(f,'units','normalized');
color=get(f,'color');

% first the data selection stuff
h_datapanel=uipanel('parent',f,...
    'title','Data Selection',...
    'tag','datapanel',...
    'fontunits','normalized',...
    'units','normalized',...
    'position',[0 0.75 0.495 0.245],...
    'backgroundcolor',color,...
    'bordertype','beveledout',...
    'borderwidth',2);
    
% second the profile criteria
h_profilepanel=uipanel('parent',f,...
    'title','Profile Criteria',...
    'tag','profilepanel',...
    'fontunits','normalized',...
    'units','normalized',...
    'position',[0 0.5 0.495 0.245],...
    'backgroundcolor',color,...
    'bordertype','beveledout',...
    'borderwidth',2);

% now the map - azimuthal projection
h_map=axesm('stereo',...
    'origin',[-a.event.centroidlat ...
    (a.event.centroidlon+180)-((a.event.centroidlon+180)>180)*360],...
    'frame','on',...
    'ffill',200,...
    'flinewidth',0.5,...
    'ffacecolor',[0.7 0.7 1]);
axis off;
set(h_map,'parent',f,...
    'tag','map',...
    'units','normalized',...
    'position',[0.5 0.5 0.5 0.5]);
% plot the land
landareas=shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
    'parent',h_map,...
    'tag','mapland',...
    'linewidth',1,...
    'facecolor',[0.7 1 0.7]);
% draw azi lines every 30deg (dotted)
plot([-2 0 -1.73205080756888 -1.73205080756888 -1  1; 
       2 0  1.73205080756888  1.73205080756888  1 -1],...
    [0 -2 -1  1 -1.73205080756888 -1.73205080756888; 
     0  2  1 -1  1.73205080756888  1.73205080756888],'k:',...
     'tag','mapazigrid',...
     'parent',h_map);
% draw dist lines every 30deg (dotted)
onedeg=(0:pi/180:2*pi)'; s1d=2*sin(onedeg); c1d=2*cos(onedeg);
plot([c1d/3 c1d/1.5],[s1d/3 s1d/1.5],'k:',...
    'tag','mapdistgrid',...
    'parent',h_map);
% make map overlay that can be selected
h_mapoverlay=axes('parent',f,...
    'tag','mapoverlay',...
    'color','none',...
    'xcolor',color,...
    'ycolor',color,...
    'xlim',[-2 2],...
    'ylim',[-2 2],...
    'xtick',[],...
    'ytick',[],...
    'xaxislocation','top',...
    'yaxislocation','right',...
    'units','normalized',...
    'position',[0.5 0.5 0.5 0.5],...
    'dataaspectratio',[1 1 1],...
    'buttondownfcn',@map_region_select);
% default region patch
azmin=option.region.azmin;
if(option.region.azmin>option.region.azmax)
    azmax=option.region.azmax+360;
else
    azmax=option.region.azmax;
end
daz=(azmax-azmin)/100;
[x,y]=da2xy([option.region.dmin*ones(1,101)...
    option.region.dmax*ones(1,101)],...
    [azmin:daz:azmax azmax:-daz:azmin]);
patch(x,y,'r',...
    'facealpha',0.5,...
    'edgecolor','r',...
    'edgealpha',0.8,...
    'tag','mappatch',...
    'parent',h_map);

% data plots
h_phasepanel=uipanel('parent',f,...
    'title','Phase',...
    'tag','phasepanel',...
    'fontunits','normalized',...
    'units','normalized',...
    'position',[0 0 0.33 0.495],...
    'backgroundcolor',color,...
    'bordertype','beveledout',...
    'borderwidth',2);
h_grouppanel=uipanel('parent',f,...
    'title','Group',...
    'tag','grouppanel',...
    'fontunits','normalized',...
    'units','normalized',...
    'position',[0.333 0 0.33 0.495],...
    'backgroundcolor',color,...
    'bordertype','beveledout',...
    'borderwidth',2);
h_amppanel=uipanel('parent',f,...
    'title','Amplitude',...
    'tag','amppanel',...
    'fontunits','normalized',...
    'units','normalized',...
    'position',[0.666 0 0.33 0.495],...
    'backgroundcolor',color,...
    'bordertype','beveledout',...
    'borderwidth',2);
% profile plots


% add figure handles, tables, defaults to guidata
data=guidata(f);

fields=fieldnames(a);
for i=1:numel(fields)
    data.(fields{i})=a.(fields{i});
end
fields=fieldnames(option);
for i=1:numel(fields)
    data.option.(fields{i})=option.(fields{i});
end
handles=guihandles(f);
fields=fieldnames(handles);
for i=1:numel(fields)
    data.(fields{i})=handles.(fields{i}); 
end
guidata(f,data);

end



function [a]=get_profiles(a)

end



function map_region_select(src,evt)
% deal with map gui input

% get data
data=guidata(src);

% get user input
[x,y,button]=ginput(2);

% convert to dist,az
az=sort(mod(atan2(x,y)*180./pi,360));
if((az(2)-az(1))>180 && all(button==1))
    data.option.region.azmin=az(2);
    data.option.region.azmax=az(1);
else
    data.option.region.azmin=az(1);
    data.option.region.azmax=az(2);
end
dist=sort(max(90,180-sqrt(x.^2+y.^2)*45));
data.option.region.dmin=dist(1);
data.option.region.dmax=dist(2);

% save data
guidata(src,data);

% update region selection
update_region_select(src,evt);

end



function update_region_select(src,evt)
% fixes bad region bounds and updates textboxes and map patch

% get data
data=guidata(src);

% fix region selection
if(data.option.region.dmin<90); data.option.region.dmin=90; end
if(data.option.region.dmax<90); data.option.region.dmax=90; end
if(data.option.region.dmin>data.option.region.dmax)
    [data.option.region.dmin,data.option.region.dmax]=...
        deal(data.option.region.dmax,data.option.region.dmin);
end
data.option.region.azmin=mod(data.option.region.azmin,360);
data.option.region.azmax=mod(data.option.region.azmax,360);

% update text boxes


% build new patch for map
azmin=data.option.region.azmin;
if(data.option.region.azmin>data.option.region.azmax)
    azmax=data.option.region.azmax+360;
else
    azmax=data.option.region.azmax;
end
daz=(azmax-azmin)/100;
[x,y]=da2xy([data.option.region.dmin*ones(1,101)...
    data.option.region.dmax*ones(1,101)],...
    [azmin:daz:azmax azmax:-daz:azmin]);

% kill old patch
delete(data.mappatch);

% apply new patch
data.mappatch=patch(x,y,'r',...
    'edgecolor','r',...
    'edgealpha',0.8,...
    'facealpha',0.5,...
    'parent',data.map);

% save data
guidata(src,data);

% update selected data
update_selected_data(src,evt);

end



function [x,y]=da2xy(d,a)
% (dist,azi) to (x,y)

x=(180-d).*sind(a)./45;
y=(180-d).*cosd(a)./45;

end



function update_selected_data(src,evt)

% get data
data=guidata(src);

% find records that meet all criteria


% save data
guidata(src,data);

% update plots
update_plots(src,evt);

end



function update_plots(src,evt)

% get data
data=guidata(src);

% update data plots


% save data
guidata(src,data);

% update profile plots
update_profile_plots(src,evt);

end



function update_profile_plots(src,evt)

% get data
data=guidata(src);

% do something with evt
if(~isempty(evt))
    error('Something went wrong here.');
end

% update plots


% save data
guidata(src,data);

end


% what subfunctions will we need
% 1  - n - spawn figure layout
% 2  - n - update subplots for profile data
% 3  - n - update subplots for selected data
% 4  - n - update data on map
% 5  - n - update patch on map
% 6  - n - update profile constraints
% 7  - n - update selected data
% 7  - n - update component list
% 8  - n - update stream list
% 9  - n - update station list
% 10 - n - update network list
% 11 - n - update region selection
% 12 - n - handle clicks on map
% 13 - n - update snr
% 14 - n - update cluster




