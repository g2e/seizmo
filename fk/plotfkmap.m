function [varargout]=plotfkmap(map,varargin)
%PLOTFKMAP    Plots the frequency-wavenumber output from FKMAP
%
%    Usage:    h=plotfkmap(map)
%              h=plotfkmap(map,fgcolor,bgcolor)
%              h=plotfkmap(map,fgcolor,bgcolor,h)
%
%    Description: H=PLOTFKMAP(MAP) plots a slowness map using the struct
%     MAP which was output from FKMAP.  See FKMAP for details on the
%     struct.  This is mainly so you can save the results and replot them
%     later (because FKMAP is quite slow).  H is the handle to the axes
%     that the map was plotted in.
%
%     H=PLOTFKMAP(MAP,FGCOLOR,BGCOLOR) specifies the foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR and
%     'k' for BGCOLOR.  Note that if one is specified and the other is not,
%     an opposing color is found using INVERTCOLOR.  The color scale is
%     also changed so the noise clip is at BGCOLOR.
%
%     H=PLOTFKMAP(MAP,FGCOLOR,BGCOLOR,H) sets the axes that the map is
%     drawn in.  This is useful for subplots, guis, etc.
%
%    Notes:
%
%    Examples:
%     Show slowness map for a dataset at about 50s periods:
%      map=fkmap(data,50,201,[1/51 1/49]);
%      plotfkmap(map);
%
%    See also: FKMAP, FKARF, UPDATEFKMAP

%     Version History:
%        May   4, 2010 - initial version
%        May  11, 2010 - updated for struct changes, got polar working,
%                        coloring options, option for setting axes handle
%        May  21, 2010 - display period rather than frequency
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  21, 2010 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% check fk struct
error(chkfkstruct(map));

% don't allow array/volume
if(~isscalar(map) || map.volume)
    error('seizmo:plotfkmap:badInput',...
        'MAP must be a scalar fk struct and not a volume!');
end

% plotting function call depends on polar
if(map.polar)
    [ax]=plotfkpolarmap(map,varargin{:});
else % cartesian
    [ax]=plotfkcartmap(map,varargin{:});
end
if(nargout); varargout{1}=ax; end

end

function ax=plotfkpolarmap(map,fgcolor,bgcolor,ax)

% check colors
if(nargin<2); fgcolor='w'; bgcolor='k'; end
if(nargin<3)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end
if(nargin<4)
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
if(nargin<4 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
end

% pertinent info
fmin=min(map.z);
fmax=max(map.z);
smax=max(abs(map.y));

% get root defaults
defaulttextcolor=get(0,'defaulttextcolor');
defaultaxescolor=get(0,'defaultaxescolor');
defaultaxesxcolor=get(0,'defaultaxesxcolor');
defaultaxesycolor=get(0,'defaultaxesycolor');
defaultaxeszcolor=get(0,'defaultaxeszcolor');
defaultpatchfacecolor=get(0,'defaultpatchfacecolor');
defaultpatchedgecolor=get(0,'defaultpatchedgecolor');
defaultlinecolor=get(0,'defaultlinecolor');
defaultsurfaceedgecolor=get(0,'defaultsurfaceedgecolor');

% set root defaults
set(0,'defaulttextcolor',fgcolor);
set(0,'defaultaxescolor',bgcolor);
set(0,'defaultaxesxcolor',fgcolor);
set(0,'defaultaxesycolor',fgcolor);
set(0,'defaultaxeszcolor',fgcolor);
set(0,'defaultpatchfacecolor',bgcolor);
set(0,'defaultpatchedgecolor',fgcolor);
set(0,'defaultlinecolor',fgcolor);
set(0,'defaultsurfaceedgecolor',fgcolor);

% initialize polar plot
ph=polar([0 2*pi],[0 smax]);
%ph=mmpolar([0 2*pi],[0 smax],'style','compass',...
%    'backgroundcolor','k','bordercolor','w');

% adjust to proper orientation
axis('ij');
delete(ph);
view([-90 90]);
hold on;

% get cartesian coords
nx=numel(map.x);
ny=numel(map.y);
[x,y]=pol2cart(pi/180*map.x(ones(ny,1),:),map.y(:,ones(nx,1)));

% plot polar grid
pcolor(x,y,map.response);

% add title color etc
hold off;
title({['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/fmin) ' to ' num2str(1/fmax) 'Sec'] ...
    '' '' ''},...
    'fontweight','bold','color',fgcolor);
set(ax,'clim',[-12 0]);
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
c=colorbar('eastoutside',...
    'fontweight','bold','xcolor',fgcolor,'ycolor',fgcolor);
set(c,'xaxislocation','top');
xlabel(c,'dB','fontweight','bold','color',fgcolor)

% reset root values
set(0,'defaulttextcolor',defaulttextcolor);
set(0,'defaultaxescolor',defaultaxescolor);
set(0,'defaultaxesxcolor',defaultaxesxcolor);
set(0,'defaultaxesycolor',defaultaxesycolor);
set(0,'defaultaxeszcolor',defaultaxeszcolor);
set(0,'defaultpatchfacecolor',defaultpatchfacecolor);
set(0,'defaultpatchedgecolor',defaultpatchedgecolor);
set(0,'defaultlinecolor',defaultlinecolor);
set(0,'defaultsurfaceedgecolor',defaultsurfaceedgecolor);

end

function ax=plotfkcartmap(map,fgcolor,bgcolor,ax)

% check colors
if(nargin<2); fgcolor='w'; bgcolor='k'; end
if(nargin<3)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end
if(nargin<4)
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
if(nargin<4 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
end

% get pertinent info
fmin=min(map.z);
fmax=max(map.z);
smax=max(max(abs(map.x)),max(abs(map.y)));

% first plot the map
imagesc(map.x,map.y,map.response);
set(ax,'xcolor',fgcolor,'ycolor',fgcolor,'ydir','normal',...
    'color',bgcolor,'fontweight','bold','clim',[-12 0]);
hold on

% phase specific bullseye
% Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
% Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
% S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
%if(smax>=37)
%    ph=[37 27.8 24.7 13.6 8.36 4.43];
%elseif(smax>=28)
%    ph=[27.8 24.7 13.6 8.36 4.43];
%elseif(smax>=25)
%    ph=[24.7 13.6 8.36 4.43];
%elseif(smax>=14)
%    ph=[13.6 8.36 4.43 2.06];
%elseif(smax>=8.5)
%    ph=[8.36 4.43 2.06];
%elseif(smax>=4.5)
%    ph=[4.43 2.06];
%else
%    ph=2.06;
%end

% regular rings (want only 3 to 5)
pot=[0.1 0.2 0.3 0.5 1 2 3 5 10 20 50 100];
rings=ceil(smax./pot);
idx=find(rings>=3 & rings<=5,1);
ph=(1:rings(idx))*pot(idx);

% plot the bull's eye
% first the radial lines
[x,y]=circle(0,12);
[x2,y2]=circle(ph(end),12);
plot([x; x2],[y; y2],'color',fgcolor,...
    'linewidth',1,'linestyle',':','tag','bullseye');
% second are the rings
for i=ph
    [x,y]=circle(i);
    plot(x,y,'color',fgcolor,'linewidth',1,'linestyle',':',...
        'tag','bullseye');
end
hold off

% finally take care of labels/coloring/etc
title({['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/fmin) ' to ' num2str(1/fmax) 'Sec']},...
    'fontweight','bold','color',fgcolor);
xlabel('East/West Slowness (s/deg)',...
    'fontweight','bold','color',fgcolor);
ylabel('North/South Slowness (s/deg)',...
    'fontweight','bold','color',fgcolor);
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
c=colorbar('eastoutside',...
    'fontweight','bold','xcolor',fgcolor,'ycolor',fgcolor);
set(c,'xaxislocation','top');
xlabel(c,'dB','fontweight','bold','color',fgcolor)
axis equal tight;

end
