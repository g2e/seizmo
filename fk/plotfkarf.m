function [varargout]=plotfkarf(arf,varargin)
%PLOTFKARF    Plots an fk array response function
%
%    Usage:    plotfkarf(arf)
%              plotfkarf(arf,dblim)
%              plotfkarf(arf,dblim,zerodb)
%              plotfkarf(arf,dblim,zerodb,fgcolor,bgcolor)
%              plotfkarf(arf,dblim,zerodb,fgcolor,bgcolor,h)
%              ax=plotfkarf(...)
%
%    Description:
%     PLOTFKARF(ARF) plots the slowness map within the struct ARF which
%     was output from FKARF.  See FKARF for details on the struct.  This is
%     mainly so you can save the results and replot them later (because
%     FKARF is slow).
%
%     PLOTFKARF(ARF,DBLIM) sets the dB limits for coloring the response
%     info.  The default is [-12 0] for the default ZERODB (see next Usage
%     form).  If ZERODB IS 'min' or 'median', the default DBLIM is [0 12].
%     DBLIM must be a real-valued 2-element vector.
%
%     PLOTFKARF(ARF,DBLIM,ZERODB) changes what 0dB corresponds to in the
%     plot.  The allowed values are 'min', 'max', 'median', & 'abs'.  The
%     default is 'max'.
%
%     PLOTFKARF(ARF,DBLIM,ZERODB,FGCOLOR,BGCOLOR) sets the foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR and
%     'k' for BGCOLOR.  Note that if one is specified and the other is not,
%     an opposing color is found using INVERTCOLOR.  The color scale is
%     also changed so the noise clip is at BGCOLOR.
%
%     PLOTFKARF(ARF,DBLIM,ZERODB,FGCOLOR,BGCOLOR,H) sets the axes that
%     the map is drawn in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTFKARF(...) returns the axis handle of the plot.
%
%    Notes:
%
%    Examples:
%     % Show a array response function for 12 plane waves:
%     arfpolar=fkarf(stla,stlo,50,201,20,[0:30:330],1/30,true);
%     plotfkarf(arfpolar);
%
%    See also: FKMAP, FKARF, PLOTFKMAP

%     Version History:
%        May  11, 2010 - initial version (outside of FKARF)
%        May  26, 2010 - update for new plotfkmap args
%        June 16, 2010 - labels call correct axes
%        July  6, 2010 - major update for new struct
%        July 18, 2010 - fix db info, commented out nyquist ring code,
%                        fixed axis output
%        Oct.  6, 2010 - truncate title if too many ARF locations
%        Oct. 10, 2010 - all plotting commands use proper ax call, tagged
%                        as 'fkarf'
%        Jan.  7, 2011 - delete commented mmpolar lines
%        Feb. 16, 2011 - fix pcolor offset in polar plots, color code fix
%        Feb. 23, 2011 - replace polar call with wedge, fix pcolor
%                        out-of-bounds pixels
%        May  25, 2011 - commented remaining out aliasing stuff
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:50 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check fk struct
error(chkfkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf))
    error('seizmo:plotfkarf:badInput',...
        'ARF must be a scalar FKARF struct!');
end

% plotting function call depends on polar
if(arf.polar)
    ax=plotfkarfpolarmap(arf,varargin{:});
else % cartesian
    ax=plotfkarfcartmap(arf,varargin{:});
end
if(nargout); varargout{1}=ax; end

end



function ax=plotfkarfpolarmap(map,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkarf:badSTYPE',...
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
    error('seizmo:plotfkarf:badDBLIM',...
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

% set zerodb & dblim in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','fkarf');

% pertinent info
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

% get nearest neighbor station distances
%[clat,clon]=arraycenter(map.stla,map.stlo);
%[e,n]=geographic2enu(clat,clon,0,map.stla,map.stlo,0);
%tri=delaunay(e,n);
%friends=[tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
%friends=unique([min(friends,[],2) max(friends,[],2)],'rows');
%dist=vincentyinv(map.stla(friends(:,1)),map.stlo(friends(:,1)),...
%                 map.stla(friends(:,2)),map.stlo(friends(:,2)));

% last plot the nyquist rings about the plane wave locations
%for i=1:map.npw
%    snyq=snyquist(min(dist),map.f(i)); % closest 2 stations
%    [x,y]=circle(snyq);
%    x=x+map.s(i)*sin(map.baz(i)*pi/180);
%    y=y+map.s(i)*cos(map.baz(i)*pi/180);
%    plot(y,x,'r:','linewidth',2,'tag','nyquist_rings');
%end

% create title
if(map.npw<=5)
    titstr=cell(map.npw,1);
    for i=1:map.npw
        titstr{i}=sprintf('SLOWNESS: %gs/deg, BAZ: %gdeg, PERIOD: %gs',...
                map.s(i),map.baz(i),1/map.f(i));
    end
else
    % elimate too much title
    titstr{1}=[num2str(map.npw) ' Locations'];
end

% add title color etc
hold(ax,'off');
title(ax,['Array Response Function @ '; titstr; ...
    ['0db = ' num2str(map.normdb) 'dB']; {''}; {''}; {''}],...
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

end



function ax=plotfkarfcartmap(map,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkarf:badSTYPE',...
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
    error('seizmo:plotfkarf:badDBLIM',...
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

% set zerodb & dblim in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','fkarf');

% get pertinent info
smax=max(max(abs(map.x)),max(abs(map.y)));

% get nearest neighbor station distances
%[clat,clon]=arraycenter(map.stla,map.stlo);
%[e,n]=geographic2enu(map.stla,map.stlo,0,clat,clon,0);
%tri=delaunay(e,n);
%friends=[tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
%friends=unique([min(friends,[],2) max(friends,[],2)],'rows');
%dist=vincentyinv(map.stla(friends(:,1)),map.stlo(friends(:,1)),...
%                 map.stla(friends(:,2)),map.stlo(friends(:,2)));

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
% next the rings
for i=ph
    [x,y]=circle(i);
    plot(ax,x,y,'color',fgcolor,'linewidth',1,'linestyle',':',...
        'tag','bullseye');
end

% last plot the nyquist rings about the plane wave locations
%for i=1:map.npw
%    snyq=snyquist(min(dist),map.f(i)); % closest 2 stations
%    [x,y]=circle(snyq);
%    x=x+map.s(i)*sin(map.baz(i)*pi/180);
%    y=y+map.s(i)*cos(map.baz(i)*pi/180);
%    plot(x,y,'r:','linewidth',2,'tag','nyquist_rings');
%end

% create title
if(map.npw<=5)
    titstr=cell(map.npw,1);
    for i=1:map.npw
        titstr{i}=sprintf('SLOWNESS: %gs/deg, BAZ: %gdeg, PERIOD: %gs',...
                map.s(i),map.baz(i),1/map.f(i));
    end
else
    % elimate too much title
    titstr{1}=[num2str(map.npw) ' Locations'];
end
hold(ax,'off');

% finally take care of labels/coloring/etc
title(ax,['Array Response Function @ '; titstr; ...
    ['0 dB = ' num2str(map.normdb) 'dB']],...
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

end
