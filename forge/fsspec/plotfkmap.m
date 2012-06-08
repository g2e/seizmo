function [varargout]=plotfkmap(map,varargin)
%PLOTFKMAP    Plots the frequency-wavenumber output from FKMAP
%
%    Usage:    h=plotfkmap(map)
%              h=plotfkmap(map,dblim)
%              h=plotfkmap(map,dblim,zerodb)
%              h=plotfkmap(map,dblim,zerodb,fgcolor,bgcolor)
%              h=plotfkmap(map,dblim,zerodb,fgcolor,bgcolor,h)
%
%    Description: H=PLOTFKMAP(MAP) plots a slowness map using the struct
%     MAP which was output from FKMAP.  See FKMAP for details on the
%     struct.  This is mainly so you can save the results and replot them
%     later (because FKMAP is quite slow).  H is the handle to the axes
%     that the map was plotted in.
%
%     H=PLOTFKMAP(MAP,DBLIM) sets the dB limits for coloring the response
%     info.  The default is [-12 0] for the default ZERODB (see next Usage
%     form).  If ZERODB IS 'min' or 'median', the default DBLIM is [0 12].
%     DBLIM must be a real-valued 2-element vector.
%
%     H=PLOTFKMAP(MAP,DBLIM,ZERODB) changes what 0dB corresponds to in the
%     plot.  The allowed values are 'min', 'max', 'median', & 'abs'.  The
%     default is 'max'.
%
%     H=PLOTFKMAP(MAP,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies foreground
%     and background colors of the plot.  The default is 'w' for FGCOLOR &
%     'k' for BGCOLOR.  Note that if one is specified and the other is not,
%     an opposing color is found using INVERTCOLOR.  The color scale is
%     also changed so the noise clip is at BGCOLOR.
%
%     H=PLOTFKMAP(MAP,DBLIM,ZERODB,FGCOLOR,BGCOLOR,H) sets the axes to draw
%     in.  This is useful for subplots, guis, etc.
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
%        May  24, 2010 - labeling the top of colorbar is broken in r2009a
%        May  26, 2010 - added dblim & zerodb args, updated docs
%        June  3, 2010 - works with new fk struct format
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  26, 2010 at 10:00 GMT

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

% rescale response
switch zerodb
    case 'min'
        map.response=map.response-min(map.response(:));
        map.normdb=map.normdb+min(map.response(:));
    case 'max'
        map.response=map.response-max(map.response(:));
        map.normdb=map.normdb+max(map.response(:));
    case 'median'
        map.response=map.response-median(map.response(:));
        map.normdb=map.normdb+median(map.response(:));
    case 'abs'
        map.response=map.response+map.normdb;
        map.normdb=0;
end

% check colors
if(nargin<4); fgcolor='w'; bgcolor='k'; end
if(nargin<5)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end
if(nargin<6)
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
if(nargin<6 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
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
    ['0 dB = ' num2str(map.normdb) 'dB'] ...
    '' '' ''},...
    'fontweight','bold','color',fgcolor);
set(ax,'clim',dblim);
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
%set(c,'xaxislocation','top');
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

% set zerodb & dblim in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata);

end



function ax=plotfkcartmap(map,cmp,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check component
if(nargin==1 || isempty(cmp)); cmp='r'; end
if(~ischar(cmp))
    error('seizmo:plotfkmap:badCMP',...
        'CMP must be ''amplitude'' ''phase'' ''real'' or ''imaginary''!');
end
isphase=false;
cmp=strmatch(cmp,{'amplitude' 'absolute' 'phase' 'real' 'rl' 'imaginary'});
switch cmp(1)
    case {1 2} % amplitude
        name='Amp Spectra';
        map.response=10*log10(abs(map.response));
    case 3 % phase
        name='Phase Spectra';
        map.response=angle(map.response);
        isphase=true;
    case {4 5} % real
        name='Real Cmp';
        map.response=10*log10(abs(real(map.response)));
    case 6 % imaginary
        name='Imag Cmp';
        map.response=10*log10(abs(imag(map.response)));
    otherwise
        error('seizmo:plotfkmap:badCMP',...
            'CMP must be ''amp'' ''phase'' ''real'' or ''imag''!');
end

% default/check scaling type
if(nargin<4 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkmap:badSTYPE',...
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
    error('seizmo:plotfkmap:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);
    
% skip rescaling if phase response
if(~isphase)
    % rescale response
    switch zerodb
        case 'min'
            map.normdb=min(map.response(:));
            map.response=map.response-map.normdb;
        case 'max'
            map.normdb=max(map.response(:));
            map.response=map.response-map.normdb;
        case 'median'
            map.normdb=median(map.response(:));
            map.response=map.response-map.normdb;
        case 'abs'
            map.normdb=0;
    end
    cblabel='dB';
    zdbtit=[name '   0 dB = ' num2str(map.normdb) 'dB'];
else
    dblim=[-pi pi];
    cblabel='radians';
    zdbtit=name;
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
if(nargin<7)
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
end

% get pertinent info
fmin=min(map.z);
fmax=max(map.z);
smax=max(max(abs(map.x)),max(abs(map.y)));

% first plot the map
imagesc(map.x,map.y,map.response);
set(ax,'xcolor',fgcolor,'ycolor',fgcolor,'ydir','normal',...
    'color',bgcolor,'fontweight','bold','clim',dblim);
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
    ['Period Range:    ' num2str(1/fmin) ' to ' num2str(1/fmax) 'Sec'] ...
    zdbtit},...
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
%set(c,'xaxislocation','top');
xlabel(c,cblabel,'fontweight','bold','color',fgcolor)
axis equal tight;

% put zerodb, dblim, cmp in userdata
% - this is for updatefkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
userdata.cmp=cmp;
set(ax,'userdata',userdata);

end
