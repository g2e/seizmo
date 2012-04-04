function [varargout]=plotfkmap(map,varargin)
%PLOTFKMAP    Plots the frequency-wavenumber output from FKMAP
%
%    Usage:    plotfkmap(map)
%              plotfkmap(map,dblim)
%              plotfkmap(map,dblim,zerodb)
%              plotfkmap(map,dblim,zerodb,fgcolor,bgcolor)
%              plotfkmap(map,dblim,zerodb,fgcolor,bgcolor,h)
%              ax=plotfkmap(...)
%
%    Description:
%     PLOTFKMAP(MAP) plots a slowness map using the struct MAP which was
%     output from FKMAP.  See FKMAP for details on the struct.  This is
%     mainly so you can save the results and replot them later (because
%     FKMAP is quite slow).
%
%     PLOTFKMAP(MAP,DBLIM) sets the dB limits for coloring the beam info.
%     The default is [-12 0] for the default ZERODB (see next Usage form).
%     If ZERODB IS 'min' or 'median', the default DBLIM is [0 12].  DBLIM
%     must be a real-valued 2-element vector.
%
%     PLOTFKMAP(MAP,DBLIM,ZERODB) changes what 0dB corresponds to in the
%     plot.  The allowed values are 'min', 'max', 'median', & 'abs'.  The
%     default is 'max'.
%
%     PLOTFKMAP(MAP,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR & 'k'
%     for BGCOLOR.  Note that if one is specified and the other is not, an
%     opposing color is found using INVERTCOLOR.  The color scale is also
%     changed so the noise clip is at BGCOLOR.
%
%     PLOTFKMAP(MAP,DBLIM,ZERODB,FGCOLOR,BGCOLOR,H) sets the axes to draw
%     in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTFKMAP(...) returns the axis handle of the plot.
%
%    Notes:
%
%    Examples:
%     % Show slowness map for a dataset at about 50s periods:
%     map=fkmap(data,50,201,[1/51 1/49]);
%     plotfkmap(map);
%
%    See also: FKMAP, PLOTFKARF, UPDATEFKMAP, FKFREQSLIDE, FKFRAMESLIDE

%     Version History:
%        May   4, 2010 - initial version
%        May  11, 2010 - updated for struct changes, got polar working,
%                        coloring options, option for setting axes handle
%        May  21, 2010 - display period rather than frequency
%        May  24, 2010 - labeling the top of colorbar is broken in r2009a
%        May  26, 2010 - added dblim & zerodb args, updated docs
%        June 16, 2010 - labels call correct axes, update see also section
%        July  6, 2010 - major update for new struct
%        Oct. 10, 2010 - all plotting functions use proper ax calls, tagged
%                        plots as 'fkmap'
%        Jan.  7, 2011 - delete commented mmpolar lines
%        Feb. 16, 2011 - fix pcolor offset in polar plots, color code fix
%        Feb. 23, 2011 - replace polar call with wedge, fix pcolor
%                        out-of-bounds pixels
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:50 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

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



function ax=plotfkpolarmap(map,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkmap:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotfkmap:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale beam
switch zerodb
    case 'min'
        map.normdb=map.normdb+min(map.beam(:));
        map.beam=map.beam-min(map.beam(:));
    case 'max'
        map.normdb=map.normdb+max(map.beam(:));
        map.beam=map.beam-max(map.beam(:));
    case 'median'
        map.normdb=map.normdb+median(map.beam(:));
        map.beam=map.beam-median(map.beam(:));
    case 'abs'
        map.beam=map.beam+map.normdb;
        map.normdb=0;
end

% check colors
if(nargin<4);
    fgcolor='w'; bgcolor='k';
elseif(nargin<5)
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
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(nargin<6 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
end

% pertinent info
fmin=min(map.freq);
fmax=max(map.freq);
smax=max(abs(map.y));

% initialize wedge plot
ph=wedge(ax,0,smax,'rcolor',fgcolor,'azcolor',fgcolor,...
                   'backgroundcolor',bgcolor,'gridcolor',fgcolor);

% clear and hold
delete(ph);
hold(ax,'on');

% offset for pcolor business
degstep=map.x(2)-map.x(1);
map.x=[map.x-degstep/2 map.x(end)+degstep/2];
slowstep=map.y(2)-map.y(1);
map.y=[map.y(1); map.y(2:end)-slowstep/2; map.y(end)]; % truncate 0 & smax

% get cartesian coords
nx=numel(map.x);
ny=numel(map.y);
[y,x]=pol2cart(pi/180*map.x(ones(ny,1),:),map.y(:,ones(nx,1)));

% plot polar grid
ph=pcolor(ax,x,y,double(map.beam([1:end end],[1:end end])));
set(ph,'clipping','off');

% add title color etc
hold(ax,'off');
title(ax,{['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/fmin) ' to ' num2str(1/fmax) 'Sec'] ...
    ['0 dB = ' num2str(map.normdb) 'dB'] ...
    '' '' ''},...
    'fontweight','bold','color',fgcolor);
set(ax,'clim',dblim);
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
c=colorbar('eastoutside','peer',ax,...
    'fontweight','bold','xcolor',fgcolor,'ycolor',fgcolor);
%set(c,'xaxislocation','top');
xlabel(c,'dB','fontweight','bold','color',fgcolor);

% set zerodb & dblim in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','fkmap');

end



function ax=plotfkcartmap(map,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkmap:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotfkmap:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale beam
switch zerodb
    case 'min'
        map.normdb=map.normdb+min(map.beam(:));
        map.beam=map.beam-min(map.beam(:));
    case 'max'
        map.normdb=map.normdb+max(map.beam(:));
        map.beam=map.beam-max(map.beam(:));
    case 'median'
        map.normdb=map.normdb+median(map.beam(:));
        map.beam=map.beam-median(map.beam(:));
    case 'abs'
        map.beam=map.beam+map.normdb;
        map.normdb=0;
end

% check colors
if(nargin<4);
    fgcolor='w'; bgcolor='k';
elseif(nargin<5)
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
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(nargin<6 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
end

% get pertinent info
fmin=min(map.freq);
fmax=max(map.freq);
smax=max(max(abs(map.x)),max(abs(map.y)));

% first plot the map
imagesc(map.x,map.y,map.beam,'parent',ax);
set(ax,'xcolor',fgcolor,'ycolor',fgcolor,'ydir','normal',...
    'color',bgcolor,'fontweight','bold','clim',dblim);
hold(ax,'on');

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
plot(ax,[x; x2],[y; y2],'color',fgcolor,...
    'linewidth',1,'linestyle',':','tag','bullseye');
% second are the rings
for i=ph
    [x,y]=circle(i);
    plot(ax,x,y,'color',fgcolor,'linewidth',1,'linestyle',':',...
        'tag','bullseye');
end
hold(ax,'off');

% finally take care of labels/coloring/etc
title(ax,{['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/fmin) ' to ' num2str(1/fmax) 'Sec'] ...
    ['0 dB = ' num2str(map.normdb) 'dB']},...
    'fontweight','bold','color',fgcolor);
xlabel(ax,'East/West Slowness (s/deg)',...
    'fontweight','bold','color',fgcolor);
ylabel(ax,'North/South Slowness (s/deg)',...
    'fontweight','bold','color',fgcolor);
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
c=colorbar('eastoutside','peer',ax,...
    'fontweight','bold','xcolor',fgcolor,'ycolor',fgcolor);
%set(c,'xaxislocation','top');
xlabel(c,'dB','fontweight','bold','color',fgcolor);
axis(ax,'equal','tight');

% set zerodb & dblim in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','fkmap');

end
