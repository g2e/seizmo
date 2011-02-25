function [varargout]=wedge(varargin)
%WEDGE    Polar coordinate wedge plot
%
%    Usage:    wedge(azimuth,radius)
%              wedge(azimuth,radius,linespec)
%              wedge(azimuth,radius,pcolormatrix)
%              wedge(...,'prop1',val1,'prop2',val2,...)
%              wedge(ax,...)
%              h=wedge(...)
%              s=wedge(ax)
%              wedge(azimuth,radius,s)
%              val=wedge(ax,'prop')
%              [val1,val2,...]=wedge(ax,{'prop1' 'prop2' ...})
%
%    Description:
%     WEDGE(AZIMUTH,RADIUS) creates a plot using azimuth coordinates in
%     degrees and radius.  Zero degrees points North (up) and increasing
%     degrees proceeds clockwise.
%
%     WEDGE(AZIMUTH,RADIUS,LINESPEC) uses the linestyle specified in string
%     LINESPEC.  See PLOT for a description of legal linestyles.
%
%     WEDGE(AZIMUTH,RADIUS,PCOLORMATRIX) passes inputs to pcolor as a way
%     to plot images in polar coordinates.  AZIMUTH & RADIUS must be
%     vectors or matrices of size NRxNAZ where NR is the number of radii &
%     NAZ is the number of azimuths.  MATRIX must be NRxNAZ in size.
%     Inputs are shifted so that the pcolor "image" pixels are at the
%     correct positions.
%
%     WEDGE(...,'PROP1',VAL1,'PROP2',VAL2,...) alters specific wedge &
%     lineseries (or surface) properties.  Available wedge properties are:
%      PROPERTY           DEFAULT
%      --------           -------
%      axeslinewidth      0.5       % linewidth of wedge axes
%      axis               'on'      % visibility of wedge axes
%      azaxislocation     'out'     % in/out 
%      azcolor            'k'       % color of azimuthal axis+ticks+labels
%      azdir              'cw'      % cw/ccw (note behavior w/ azdir)
%      azgrid             'on'      % visibility of azimuthal grid lines
%      azlim              []        % defaults to [0 360] or [0 2*pi]
%      azoffset           0         % 0 is vert, affected by azdir, azunits
%      aztick             []        % suitable values found by default
%      aztickdelta        []        % suitable delta found by default
%      azticklabel        []        % uses aztick values
%      azticklabelhalign  'center'  % center/left/right/auto
%      azticklabeloffset  20        % radial offset in pixels
%      azticklabelvalign  'middle'  % middle/top/cap/baseline/bottom
%      azticklabelvisible 'on'      % visibility of azimuthal tick labels
%      aztickvisible      'on'      % visibility of azimuthal ticks
%      azunits            'degrees' % degrees/radians
%      backgroundcolor    'w'       % color of wedge plot background patch
%      box                'on'      % boxed wedge plot
%      fontangle          ----      % uses get(0,'defaulttextfontangle')
%      fontname           ----      % uses get(0,'defaulttextfontname')
%      fontsize           ----      % uses get(0,'defaulttextfontsize')
%      fontunits          ----      % uses get(0,'defaulttextfontunits')
%      fontweight         ----      % uses get(0,'defaulttextfontweight')
%      gridcolor          'k'       % color of grid lines
%      gridlinestyle      ':'       % linestyle of grid lines
%      gridlinewidth      0.5       % linewidth of grid lines
%      npts               51        % number of points in azimuthal axis
%      placement          'center'  % center/origin
%      plotbox            'off'     % on/off (clip wedge during pan/zoom)
%      raxislocation      'cw'      % cw/ccw (only drawn for wedge)
%      rcolor             'k'       % color of radial axis+ticks+labels
%      rgrid              'on'      % visibility of radial grid lines
%      rlim               []        % defaults to [0 max(radius)]
%      rtick              []        % suitable values found by default
%      rtickdelta         []        % suitable values found by default
%      rticklabel         []        % uses rtick values
%      rticklabelcircaz   0         % r axis label azimuth when circle
%      rticklabelhalign   'center'  % center/left/right/auto
%      rticklabeloffset   20        % azimuthal offset in pixels
%      rticklabelvalign   'middle'  % middle/top/cap/baseline/bottom
%      rticklabelvisible  'on'      % visibility of radial tick labels
%      rticklabelunits    ''        % units string to append to labels
%      rtickvisible       'on'      % visibility of radial ticks
%      tag                ''        % tag for wedge axes
%      tickdir            'in'      % direction of ticks
%      ticklength         0.02      % fraction of rlim difference
%      z                  0         % z value for wedge axes
%     All other property/value pairs are passed on to PLOT (or PCOLOR) when
%     drawing the input data.  
%
%     WEDGE(AX,...) plots in AX instead of GCA.
%
%     H=WEDGE(...) returns the handle to the plotted object as H.
%
%     S=WEDGE(AX) or S=WEDGE returns a struct containing the properties of
%     the wedge axes specified by AX.  If no axes is given then the current
%     axes properties are output (if it is a wedge plot).
%
%     WEDGE(AZIMUTH,RADIUS,S) allows setting properties of the wedge axes
%     using the input structure S.
%
%     VAL=WEDGE(AX,'PROP') returns the value of the specified wedge
%     property.  See above for valid properties.
%
%     [VAL1,VAL2,...]=WEDGE(AX,{'PROP1' 'PROP2' ...}) returns multiply
%     properties of the wedge axes.
%
%    Notes:
%
%    Examples:
%     % Polar plot:
%     wedge(360*rand(10,1),rand(10,1),'azdir','ccw','azoffset',-90);
%
%     % Radioactive wedge:
%     wedge(axes,[],[],'azlim',[-15 15],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%     wedge(axes,[],[],'azlim',[105 135],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%     wedge(axes,[],[],'azlim',[225 255],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%
%    See also: POLAR, ROSE, COMPASS

%     Version History:
%        Feb. 22, 2011 - initial version
%        Feb. 23, 2011 - fix missing handle output, add padding for text,
%                        changed halign & offset defaults
%        Feb. 24, 2011 - pcolor option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 24, 2011 at 18:35 GMT

% todo
% - drawing
%   - minor grid/ticks
%   - clipping
% - updating wedges
% - rlabel/azlabel (separate functions)
% - multiple wedges in the same axes
% - dms stuff
% - colormap background

% find axes input
[ax,varargin,nargs]=axescheck(varargin{:});

% find current axes if none given
% - this also returns current if given was dead
if(isempty(ax)); ax=get(get(0,'currentfigure'),'currentaxes'); end

% act based on input style
if(~nargs) % struct output
    require_wedge(ax);
    varargout{1}=getappdata(ax);
elseif(nargs==1) % get property or set properties using struct
    require_wedge(ax);
    if(isstruct(varargin{1})) % set
        wedge_update(ax,wedge_check(struct2pv(varargin{1})));
    elseif(ischar(varargin{1}) || iscellstr(varargin{1})) % get
        varargin{1}=cellstr(varargin{1});
        for i=1:numel(varargin{1})
            varargout{i}=getappdata(ax,varargin{1}{i});
        end
    else
        error('seizmo:wedge:badInput',...
            'Unknown usage form (or bad wedge handle?)!');
    end
else % p/v pairs or data
    if(isnumeric(varargin{1}))
        % check inputs
        if(~isnumeric(varargin{2})) % maybe stray wedge handle + something
            error('seizmo:wedge:badInput',...
                'Unknown usage form (or bad wedge handle?)!');
        elseif(nargs>=3 && isnumeric(varargin{3}) ...
                && (~isvector(varargin{1}) ...
                || ~isvector(varargin{2}) ...
                || ~isequal(size(varargin{3}),...
                [numel(varargin{2}) numel(varargin{1})])) ...
                && ~isequal(size(varargin{1}),size(varargin{2},...
                size(varargin{3}))))
            error('seizmo:wedge:badInput',...
                'Improper wedge image call (or bad wedge handle?)!');
        end
        w=wedge_defaults();
        if(any(nargs==[3 4]) && isstruct(varargin{end}))
            varargin=[varargin(1:end-1) struct2pv(varargin{end})];
        end
        if(nargs>=3 && isnumeric(varargin{3})); lni=3; else lni=2; end
        [wpv,lpv]=wedge_splitpv(w,varargin(lni+1:end));
        if(isempty(ax) || ~ishandle(ax) || ~ishold(ax) ...
                || ~strcmp(get(ax,'createfcn'),'wedge')) % new plot
            w=appendpv2struct(w,wedge_check(wpv));
            ax=wedge_create(ax,varargin{2},w);
            if(lni==2) % lines
                h=wedge_plot(ax,varargin{1:2},lpv,w,false);
            else % image
                h=wedge_pcolor(ax,varargin{1:3},lpv,w,false);
            end
            if(nargout); varargout{1}=h; end
        elseif(isempty(wpv)) % old plot, more data
            if(lni==2) % lines
                h=wedge_plot(ax,varargin{1:2},lpv,getappdata(ax),true);
            else % image
                h=wedge_pcolor(ax,varargin{1:3},lpv,getappdata(ax),true);
            end
            if(nargout); varargout{1}=h; end
        else % update old plot, more data
            wedge_update(ax,wedge_check(wpv));
            if(lni==2) % lines
                h=wedge_plot(ax,varargin{1:2},lpv,getappdata(ax),true);
            else % image
                h=wedge_pcolor(ax,varargin{1:3},lpv,getappdata(ax),true);
            end
            if(nargout); varargout{1}=h; end
        end
    elseif(ischar(varargin{1})) % p/v
        require_wedge(ax);
        wedge_update(ax,wedge_check(varargin{1}));
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h]=wedge_plot(ax,x,y,lpv,w,held)
% convert to cartesian
% - need to account for azdir, azoffset, azunits
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
if(strcmp(w.azunits,'radians')); wsin=@sin; wcos=@cos;
else wsin=@sind; wcos=@cosd;
end
[x,y]=deal(y.*wsin(azdir*(x+w.azoffset)),y.*wcos(azdir*(x+w.azoffset)));

% plot data on top
h=plot(x,y,lpv{:},'parent',ax,'clipping',w.plotbox);

% cleanup
set(ax,'createfcn','wedge');
if(~held)
    hold(ax,'off');
    axis(ax,'off');
    set(get(ax,'xlabel'),'visible','on');
    set(get(ax,'ylabel'),'visible','on');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h]=wedge_pcolor(ax,x,y,z,lpv,w,held)
% vector 2 matrices
if(isvector(x) && isvector(y))
    nx=numel(x); ny=numel(y);
    degstep=x(2)-x(1);
    radstep=y(2)-y(1);
    x=x(:)';
    y=y(:);
    x=x(ones(ny+1,1),:);
    y=y(:,ones(1,nx+1));
else % check matrices
    naz=size(x,2);
    nr=size(y,1);
    if(numel(unique(x(:)))~=naz ...
            || ~isequal(unique(x(:)),unique(x(1,:))'))
        error('seizmo:wedge:badInput',...
            'AZIMUTH matrix not formatted properly!');
    elseif(numel(unique(y(:)))~=nr ...
            || ~isequal(unique(y(:)),unique(y(:,1))))
        error('seizmo:wedge:badInput',...
            'RADIUS matrix not formatted properly!');
    elseif(~isequal(size(x),size(y),size(z),[nr naz]))
        error('seizmo:wedge:badInput',...
            'Matrices for PCOLOR input not formatted properly!');
    end
    x=x([1:end end],:);
    y=y(:,[1:end end]);
end

% image to pcolor
x=[x-degstep/2 x(:,end)+degstep/2];
y=[y-radstep/2; y(end,:)+radstep/2];
z=z([1:end end],[1:end end]);

% clipping (sorta)
if(strcmpi(w.clipping,'on'))
    x(x<w.azlim(1))=w.azlim(1);
    x(x>w.azlim(2))=w.azlim(2);
    y(y<w.rlim(1))=w.rlim(1);
    y(y>w.rlim(2))=w.rlim(2);
end

% polar 2 cartesian
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
if(strcmp(w.azunits,'radians')); wsin=@sin; wcos=@cos;
else wsin=@sind; wcos=@cosd;
end
[x,y]=deal(y.*wsin(azdir*(x+w.azoffset)),y.*wcos(azdir*(x+w.azoffset)));

% pcolor
h=pcolor(ax,x,y,z);
set(h,'clipping',w.plotbox,lpv{:});
shading(ax,'flat');

% cleanup
set(ax,'createfcn','wedge');
if(~held)
    hold(ax,'off');
    axis(ax,'off');
    set(get(ax,'xlabel'),'visible','on');
    set(get(ax,'ylabel'),'visible','on');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax]=wedge_create(ax,y,w)
%WEDGE_CREATE    Creates wedge plot

% clear plot & set hold to on
if(ishandle(ax)); hold(ax,'off'); end
ax=newplot(ax);
hold(ax,'on');

% first determine auto rlim & tlim
if(isempty(w.rlim))
    rmax=max(abs(y(~isinf(y))));
    if(isempty(rmax)); rmax=1; end
    w.rlim=[0 rmax];
end
if(isempty(w.azlim))
    if(strcmp(w.azunits,'degrees')); w.azlim=[0 360];
    else w.azlim=[0 2*pi];
    end
else
    % force over circled to just a circle
    if(strcmp(w.azunits,'degrees'))
        if(diff(w.azlim)>360); w.azlim=w.azlim(1)+[0 360]; end
    else
        if(diff(w.azlim)>2*pi); w.azlim=w.azlim(1)+[0 2*pi]; end
    end
end

% draw axes
wedge_axes(ax,w);

% get & draw ticks
w=wedge_ticks(ax,w);

% draw tick labels
wedge_ticklabels(ax,w);

% draw grid
wedge_grid(ax,w);

% draw minorticks
wedge_minorticks(ax,w);

% draw minorgrid
wedge_minorgrid(ax,w);

% save options to appdata
fields=fieldnames(w);
for i=1:numel(fields); setappdata(ax,fields{i},w.(fields{i})); end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w]=wedge_defaults()
w.axeslinewidth=0.5;
w.axis='on';
w.azaxislocation='out'; % in/out
w.azcolor='k';
w.azdir='cw'; % cw/ccw
w.azgrid='on';
w.azlim=[];
%w.azminorgrid='off';
%w.azminortick='off';
w.azoffset=0; % 0 is vertical, affected by azdir, azunits
w.aztick=[];
w.aztickdelta=[];
w.azticklabel=[];
w.azticklabelhalign='center'; % center/left/right/auto
w.azticklabeloffset=20; % outward offset in pixels
w.azticklabelvalign='middle'; % middle/top/cap/baseline/bottom
w.azticklabelvisible='on';
w.aztickvisible='on';
w.azunits='degrees'; % degrees/radians
w.backgroundcolor='w';
w.box='on';
w.clipping='on'; % only works (sorta) for image
w.fontangle=get(0,'defaulttextfontangle');
w.fontname=get(0,'defaulttextfontname');
w.fontsize=get(0,'defaulttextfontsize');
w.fontunits=get(0,'defaulttextfontunits');
w.fontweight=get(0,'defaulttextfontweight');
w.gridcolor='k';
w.gridlinestyle=':';
w.gridlinewidth=0.5;
%w.minorgridcolor='k';
%w.minorgridlinestyle=':';
%w.minorgridlinewidth=0.5;
w.npts=51;
w.placement='center'; % center/origin
w.plotbox='off'; % on/off (clip wedge during pan/zoom)
w.raxislocation='cw'; % cw/ccw for wedge
w.rcolor='k';
w.rgrid='on';
w.rlim=[];
%w.rminorgrid='off';
%w.rminortick='off';
w.rtick=[];
w.rtickdelta=[];
w.rticklabel=[];
w.rticklabelcircaz=0; % r axis label azimuth when circle
w.rticklabelhalign='center'; % center/left/right/auto
w.rticklabeloffset=20; % outward offset in pixels
w.rticklabelvalign='middle'; % middle/top/cap/baseline/bottom
w.rticklabelvisible='on';
w.rticklabelunits='';
w.rtickvisible='on';
w.tag='';
w.tickdir='in';
w.ticklength=0.02; % normalized to rlim difference
w.z=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_axes(ax,w)

% setup for azimuth units & direction
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
cw=(3-azdir)/2; ccw=(3+azdir)/2;
if(strcmp(w.azunits,'radians')); azdeg=false; wsin=@sin; wcos=@cos;
else azdeg=true; wsin=@sind; wcos=@cosd;
end

% full circle?
if(azdeg); circ=diff(w.azlim)>=360;
else circ=diff(w.azlim)>=2*pi;
end

% hollow?
if(w.rlim(1)<=0); hollow=false; else hollow=true; end

% axes lines
axesinx=[];  axesiny=[];
axescwx=[];  axescwy=[];
axesccwx=[]; axesccwy=[];
axesoutx=wsin(azdir*(linspace(w.azlim(1),w.azlim(2),w.npts)'+w.azoffset));
axesouty=wcos(azdir*(linspace(w.azlim(1),w.azlim(2),w.npts)'+w.azoffset));
if(hollow)
    % in
    axesinx=w.rlim(1)*axesoutx;
    axesiny=w.rlim(1)*axesouty;
end
axesoutx=w.rlim(2)*axesoutx;
axesouty=w.rlim(2)*axesouty;
if(~circ)
    axescwx=[w.rlim(1)*wsin(azdir*(w.azlim(cw)+w.azoffset));...
             w.rlim(2)*wsin(azdir*(w.azlim(cw)+w.azoffset))];
    axescwy=[w.rlim(1)*wcos(azdir*(w.azlim(cw)+w.azoffset));...
             w.rlim(2)*wcos(azdir*(w.azlim(cw)+w.azoffset))];
    axesccwx=[w.rlim(1)*wsin(azdir*(w.azlim(ccw)+w.azoffset));...
              w.rlim(2)*wsin(azdir*(w.azlim(ccw)+w.azoffset))];
    axesccwy=[w.rlim(1)*wcos(azdir*(w.azlim(ccw)+w.azoffset));...
              w.rlim(2)*wcos(azdir*(w.azlim(ccw)+w.azoffset))];
end

% axes background
patch('xdata',[axesinx; axescwx; flipud(axesoutx); flipud(axesccwx)],...
    'ydata',[axesiny; axescwy; flipud(axesouty); flipud(axesccwy)],...
    'zdata',w.z(ones(numel([axesinx; axescwx; axesoutx; axesccwx]),1)),...
    'edgecolor','none','facecolor',w.backgroundcolor,'parent',ax,...
    'handlevisibility','off','tag','wedge_bgpatch','clipping',w.plotbox,...
    'visible',w.axis);

% axes
if(~isempty(axesinx))
    switch w.azaxislocation
        case 'in'
            onoff=onoffstate(w.axis);
        case 'out'
            onoff=onoffstate(w.axis,w.box);
    end
    line(axesinx,axesiny,w.z(ones(numel(axesinx),1)),...
        'parent',ax,'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_in',...
        'clipping',w.plotbox,'visible',onoff);
end
switch w.azaxislocation
    case 'in'
        onoff=onoffstate(w.axis,w.box);
    case 'out'
        onoff=onoffstate(w.axis);
end
line(axesoutx,axesouty,w.z(ones(numel(axesoutx),1)),...
    'parent',ax,'color',w.azcolor,'linewidth',w.axeslinewidth,...
    'handlevisibility','off','tag','wedge_axis_out',...
    'clipping',w.plotbox,'visible',onoff);
if(~isempty(axescwx))
    switch w.raxislocation
        case 'cw'
            onoff=onoffstate(w.axis);
        case 'ccw'
            onoff=onoffstate(w.axis,w.box);
    end
    line(axescwx,axescwy,w.z(ones(numel(axescwx),1)),...
        'parent',ax,'color',w.rcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_cw',...
        'clipping',w.plotbox,'visible',onoff);
end
if(~isempty(axesccwx))
    switch w.raxislocation
        case 'cw'
            onoff=onoffstate(w.axis,w.box);
        case 'ccw'
            onoff=onoffstate(w.axis);
    end
    line(axesccwx,axesccwy,w.z(ones(numel(axesccwx),1)),...
        'parent',ax,'color',w.rcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_ccw',...
        'clipping',w.plotbox,'visible',onoff);
end

% axes placement + 15% padding
axis(ax,'equal');
if(strcmp(w.placement,'origin'))
    axis(ax,w.rlim(2)*[-1.15 1.15 -1.15 1.15]);
else
    axis(ax,'tight');
    xrng=get(ax,'xlim');
    yrng=get(ax,'ylim');
    set(ax,'xlim',xrng+.15*[-1 1]*diff(xrng)/2,...
           'ylim',yrng+.15*[-1 1]*diff(yrng)/2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w]=wedge_ticks(ax,w)

% setup for azimuth units & directions
tdir=1; if(strcmp(w.tickdir,'out')); tdir=-1; end
azaxpo=w.rlim(1); if(strcmp(w.azaxislocation,'out')); azaxpo=w.rlim(2); end
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
cw=(3-azdir)/2; ccw=(3+azdir)/2;
if(strcmp(w.azunits,'radians')); azdeg=false; wsin=@sin; wcos=@cos;
else azdeg=true; wsin=@sind; wcos=@cosd;
end

% full circle?
if(azdeg); circ=diff(w.azlim)>=360;
else circ=diff(w.azlim)>=2*pi;
end

% hollow?
if(w.rlim(1)<=0); hollow=false; else hollow=true; end

% auto tick setup
if(isempty(w.rtick) || (isempty(w.aztick) && isempty(w.aztickdelta)))
    % radial, azimuthal length
    xlim=get(ax,'xlim');
    units=get(ax,'units');
    set(ax,'units','normalized');
    pos=get(ax,'position');
    set(ax,'units',units);
    rwidth=diff(w.rlim)*pos(3)/diff(xlim);
    if(azdeg); azwidth=diff(w.azlim)*pi/180*azaxpo*pos(3)/diff(xlim);
    else azwidth=diff(w.azlim)*azaxpo*pos(3)/diff(xlim);
    end
    
    % font size
    set(ax,'fontname',w.fontname,'fontunits',w.fontunits,...
        'fontsize',w.fontsize,'fontweight',w.fontweight,...
        'fontangle',w.fontangle);
    set(ax,'fontunits','normalized');
    fsize=get(ax,'fontsize');
    set(ax,'fontunits',w.fontunits);
end

% get auto ticks
if(isempty(w.aztick) && isempty(w.aztickdelta))
    % full auto
    if(azdeg)
        if(circ)
            [w.aztick,w.aztickdelta]=magicticks360(w.azlim,azwidth,fsize);
        else
            [w.aztick,w.aztickdelta]=magicticks(w.azlim,azwidth,fsize);
        end
    else
        [w.aztick,w.aztickdelta]=magicticks(w.azlim./pi,azwidth,fsize);
        w.aztick=w.aztick*pi;
        w.aztickdelta=w.aztickdelta*pi;
    end
    if(circ)
        switch w.azunits
            case 'degrees'
                if(w.aztick(1)+360==w.aztick(end)); w.aztick(end)=[]; end
            case 'radians'
                if(w.aztick(1)+2*pi==w.aztick(end)); w.aztick(end)=[]; end
        end
    end
elseif(isempty(w.aztick))
    % delta given
    w.aztick=w.aztickdelta...
        *(ceil(w.azlim(1)/w.aztickdelta):floor(w.azlim(2)/w.aztickdelta));
    if(circ)
        switch w.azunits
            case 'degrees'
                if(w.aztick(1)+360==w.aztick(end)); w.aztick(end)=[]; end
            case 'radians'
                if(w.aztick(1)+2*pi==w.aztick(end)); w.aztick(end)=[]; end
        end
    end
end
if(isempty(w.rtick) && isempty(w.rtickdelta))
    [w.rtick,w.rtickdelta]=magicticks(w.rlim,rwidth,fsize);
    if(~w.rtick(1)); w.rtick(1)=[]; end
elseif(isempty(w.rtick))
    % delta given
    w.rtick=w.rtickdelta*(...
        ceil(w.rlim(1)/w.rtickdelta):floor(w.rlim(2)/w.rtickdelta));
    if(~w.rtick(1)); w.rtick(1)=[]; end
end

% tick mark length
abstlen=diff(w.rlim)*w.ticklength;

% tick lines
tickinx=[];  tickiny=[];
tickcwx=[];  tickcwy=[];
tickccwx=[]; tickccwy=[];
tickoutr=w.rlim(2)+[zeros(1,numel(w.aztick)); 
                    -tdir*abstlen(ones(1,numel(w.aztick)))];
tickoutx=tickoutr.*wsin(azdir*(w.aztick([1 1],:)+w.azoffset));
tickouty=tickoutr.*wcos(azdir*(w.aztick([1 1],:)+w.azoffset));
if(hollow)
    % in
    tickinr=w.rlim(1)+[zeros(1,numel(w.aztick)); 
                       tdir*abstlen(ones(1,numel(w.aztick)))];
    tickinx=tickinr.*wsin(azdir*(w.aztick([1 1],:)+w.azoffset));
    tickiny=tickinr.*wcos(azdir*(w.aztick([1 1],:)+w.azoffset));
end
if(~circ)
    % cw/ccw
    tickcwx=w.rtick*wsin(azdir*(w.azlim(cw)+w.azoffset));
    tickcwy=w.rtick*wcos(azdir*(w.azlim(cw)+w.azoffset));
    tickcwx=[tickcwx; tickcwx+tdir*abstlen...
        *wcos(azdir*(w.azlim(cw)+w.azoffset))];
    tickcwy=[tickcwy; tickcwy-tdir*abstlen...
        *wsin(azdir*(w.azlim(cw)+w.azoffset))];
    tickccwx=w.rtick*wsin(azdir*(w.azlim(ccw)+w.azoffset));
    tickccwy=w.rtick*wcos(azdir*(w.azlim(ccw)+w.azoffset));
    tickccwx=[tickccwx; tickccwx-tdir*abstlen...
        *wcos(azdir*(w.azlim(ccw)+w.azoffset))];
    tickccwy=[tickccwy; tickccwy+tdir*abstlen...
        *wsin(azdir*(w.azlim(ccw)+w.azoffset))];
end

% draw ticks
switch w.azaxislocation
    case 'in'
        onoff=onoffstate(w.axis,w.box,w.aztickvisible);
    case 'out'
        onoff=onoffstate(w.axis,w.aztickvisible);
end
line(tickoutx,tickouty,w.z(ones(2,numel(w.aztick))),'parent',ax,...
    'color',w.azcolor,'linewidth',w.axeslinewidth,...
    'handlevisibility','off','tag','wedge_tick_out',...
    'clipping',w.plotbox,'visible',onoff);
if(~isempty(tickinx))
    switch w.azaxislocation
        case 'in'
            onoff=onoffstate(w.axis,w.aztickvisible);
        case 'out'
            onoff=onoffstate(w.axis,w.box,w.aztickvisible);
    end
    line(tickinx,tickiny,w.z(ones(2,numel(w.aztick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_in',...
        'clipping',w.plotbox,'visible',onoff);
end
if(~isempty(tickcwx))
    switch w.raxislocation
        case 'cw'
            onoff=onoffstate(w.axis,w.rtickvisible);
        case 'ccw'
            onoff=onoffstate(w.axis,w.box,w.rtickvisible);
    end
    line(tickcwx,tickcwy,w.z(ones(2,numel(w.rtick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_cw',...
        'clipping',w.plotbox,'visible',onoff);
end
if(~isempty(tickccwx))
    switch w.raxislocation
        case 'cw'
            onoff=onoffstate(w.axis,w.box,w.rtickvisible);
        case 'ccw'
            onoff=onoffstate(w.axis,w.rtickvisible);
    end
    line(tickccwx,tickccwy,w.z(ones(2,numel(w.rtick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_ccw',...
        'clipping',w.plotbox,'visible',onoff);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_grid(ax,w)

% setup for azimuth units & direction
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
if(strcmp(w.azunits,'radians')); wsin=@sin; wcos=@cos;
else wsin=@sind; wcos=@cosd;
end

% radial grid (circles)
%x=r.*wsin(azdir*(az+w.azoffset))
%y=r.*wcos(azdir*(az+w.azoffset))
azpnts=linspace(w.azlim(1),w.azlim(2),w.npts)';
rgridx=w.rtick(ones(w.npts,1),:) ...
    .*wsin(azdir*(azpnts(:,ones(1,numel(w.rtick)))+w.azoffset));
rgridy=w.rtick(ones(w.npts,1),:) ...
    .*wcos(azdir*(azpnts(:,ones(1,numel(w.rtick)))+w.azoffset));

% azimuthal grid (spokes)
rpnts=w.rlim';
azgridx=rpnts(:,ones(1,numel(w.aztick))) ...
    .*wsin(azdir*(w.aztick(ones(2,1),:)+w.azoffset));
azgridy=rpnts(:,ones(1,numel(w.aztick))) ...
    .*wcos(azdir*(w.aztick(ones(2,1),:)+w.azoffset));

% plot radial grid
line(rgridx,rgridy,w.z(ones(w.npts,numel(w.rtick))),'parent',ax,...
    'color',w.gridcolor,'linewidth',w.gridlinewidth,...
    'linestyle',w.gridlinestyle,'handlevisibility','off',...
    'visible',onoffstate(w.axis,w.rgrid),...
    'tag','wedge_grid_r','clipping',w.plotbox);

% plot azimuthal grid
line(azgridx,azgridy,w.z(ones(2,numel(w.aztick))),'parent',ax,...
    'color',w.gridcolor,'linewidth',w.gridlinewidth,...
    'linestyle',w.gridlinestyle,'handlevisibility','off',...
    'visible',onoffstate(w.axis,w.azgrid),...
    'tag','wedge_grid_az','clipping',w.plotbox);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_ticklabels(ax,w)

% setup for azimuth units & directions
azdir=1; if(strcmp(w.azdir,'ccw')); azdir=-1; end
cw=(3-azdir)/2; ccw=(3+azdir)/2;
if(strcmp(w.azunits,'radians')); azdeg=false; wsin=@sin; wcos=@cos;
else azdeg=true; wsin=@sind; wcos=@cosd;
end

% label units
if(azdeg)
    azticklabelunits='^o';
else
    azticklabelunits='\pi';
end

% full circle?
if(azdeg); circ=diff(w.azlim)>=360;
else circ=diff(w.azlim)>=2*pi;
end

% hollow?
if(w.rlim(1)<=0); hollow=false; else hollow=true; end

% auto labels
if(isempty(w.azticklabel))
    if(azdeg)
        w.azticklabel=cellstr(num2str(w.aztick'));
    else
        w.azticklabel=cellstr(num2str(w.aztick'./pi));
    end
end
if(isempty(w.rticklabel))
    w.rticklabel=cellstr(num2str(w.rtick'));
end

% trim spaces
w.azticklabel=strtrim(w.azticklabel);
w.rticklabel=strtrim(w.rticklabel);

% radial tick label positioning
if(circ)
    % labels are at position specified by rticklabelcircaz
    rlabelx=w.rtick.*wsin(azdir*(w.rticklabelcircaz+w.azoffset));
    rlabely=w.rtick.*wcos(azdir*(w.rticklabelcircaz+w.azoffset));
else
    switch w.raxislocation
        case 'cw'
            rlabelx=w.rtick.*wsin(azdir*(w.azlim(cw)+w.azoffset));
            rlabely=w.rtick.*wcos(azdir*(w.azlim(cw)+w.azoffset));
        case 'ccw'
            rlabelx=w.rtick.*wsin(azdir*(w.azlim(ccw)+w.azoffset));
            rlabely=w.rtick.*wcos(azdir*(w.azlim(ccw)+w.azoffset));
    end
end

% azimuthal tick label positioning
switch w.azaxislocation
    case 'in'
        if(hollow)
            azlabelx=w.rlim(1)*wsin(azdir*(w.aztick+w.azoffset));
            azlabely=w.rlim(1)*wcos(azdir*(w.aztick+w.azoffset));
        else % do not plot labels
            azlabelx=[];
            azlabely=[];
        end
    case 'out'
        azlabelx=w.rlim(2)*wsin(azdir*(w.aztick+w.azoffset));
        azlabely=w.rlim(2)*wcos(azdir*(w.aztick+w.azoffset));
end

% auto horizontal alignment
rhalign=cellstr(w.rticklabelhalign);
rhalign=rhalign(ones(numel(w.rtick),1));
if(strcmp(w.rticklabelhalign,'auto'))
    % auto pushes text outwards from plot
    % cw &  +y => right
    %        0 => center
    %       -y => left
    % ccw & +y => left
    %        0 => center
    %       -y => right
    for i=1:numel(rlabelx)
        switch w.raxislocation
            case 'cw'
                if(rlabely(i)>0)
                    rhalign{i}='right';
                elseif(rlabely(i)<0)
                    rhalign{i}='left';
                else
                    rhalign{i}='center';
                end
            case 'ccw'
                if(rlabely(i)>0)
                    rhalign{i}='left';
                elseif(rlabely(i)<0)
                    rhalign{i}='right';
                else
                    rhalign{i}='center';
                end
        end
    end
end
azhalign=cellstr(w.azticklabelhalign);
azhalign=azhalign(ones(numel(w.aztick),1));
if(strcmp(w.azticklabelhalign,'auto') && ~isempty(azlabelx))
    % auto pushes text outwards from plot
    % out & +x => left
    %        0 => center
    %       -x => right
    % in &  +x => right
    %        0 => center
    %       -x => left
    for i=1:numel(azlabelx)
        switch w.azaxislocation
            case 'in'
                if(azlabelx(i)>eps('single'))
                    azhalign{i}='right';
                elseif(azlabelx(i)<-eps('single'))
                    azhalign{i}='left';
                else
                    azhalign{i}='center';
                end
            case 'out'
                if(azlabelx(i)>eps('single'))
                    azhalign{i}='left';
                elseif(azlabelx(i)<-eps('single'))
                    azhalign{i}='right';
                else
                    azhalign{i}='center';
                end
        end
    end
end

% plot labels
rh=nan(numel(w.rtick),1);
for i=1:numel(w.rtick)
    w.rticklabel{i}=[w.rticklabel{i} w.rticklabelunits];
    rh(i)=text(rlabelx(i),rlabely(i),w.z,w.rticklabel{i},'parent',ax,...
        'horizontalalignment',rhalign{i},'clipping',w.plotbox,...
        'visible',onoffstate(w.axis,w.rticklabelvisible),...
        'verticalalignment',w.rticklabelvalign,'color',w.rcolor,...
        'tag','wedge_rlabels','fontunits',w.fontunits,...
        'fontsize',w.fontsize,'fontname',w.fontname,...
        'fontangle',w.fontangle,'fontweight',w.fontweight);
end
if(~isempty(azlabelx))
    azh=nan(numel(w.aztick),1);
    for i=1:numel(w.aztick)
        if(~azdeg && strcmp(w.azticklabel{i},'0'))
            % no units added on
        elseif(~azdeg && strcmp(w.azticklabel{i},'1'))
            % drop the 1
            w.azticklabel{i}=azticklabelunits;
        else
            w.azticklabel{i}=[w.azticklabel{i} azticklabelunits];
        end
        azh(i)=text(azlabelx(i),azlabely(i),w.z,w.azticklabel{i},...
            'parent',ax,'horizontalalignment',azhalign{i},...
            'clipping',w.plotbox,...
            'visible',onoffstate(w.axis,w.azticklabelvisible),...
            'verticalalignment',w.azticklabelvalign,'color',w.azcolor,...
            'tag','wedge_azlabels','fontunits',w.fontunits,...
            'fontsize',w.fontsize,'fontname',w.fontname,...
            'fontangle',w.fontangle,'fontweight',w.fontweight);
    end
end

% offset labels
if(~circ) % offset radial tick labels only if not circle
    for i=1:numel(w.rtick)
        set(rh(i),'units','pixels');
        pos=get(rh(i),'position');
        % modify position
        switch w.raxislocation
            case 'cw'
                pos(1)=pos(1)-w.rticklabeloffset...
                    *wcos(azdir*(w.azlim(cw)+w.azoffset));
                pos(2)=pos(2)+w.rticklabeloffset...
                    *wsin(azdir*(w.azlim(cw)+w.azoffset));
            case 'ccw'
                pos(1)=pos(1)+w.rticklabeloffset...
                    *wcos(azdir*(w.azlim(ccw)+w.azoffset));
                pos(2)=pos(2)-w.rticklabeloffset...
                    *wsin(azdir*(w.azlim(ccw)+w.azoffset));
        end
        set(rh(i),'position',pos);
        set(rh(i),'units','data');
        % fix matlab bug
        pos=get(rh(i),'position');
        pos(3)=w.z;
        set(rh(i),'position',pos);
    end
end
if(~isempty(azlabelx))
    for i=1:numel(w.aztick)
        set(azh(i),'units','pixels');
        pos=get(azh(i),'position');
        % modify position
        switch w.azaxislocation
            case 'in'
                pos(1)=pos(1)-w.azticklabeloffset...
                    *wsin(azdir*(w.aztick(i)+w.azoffset));
                pos(2)=pos(2)-w.azticklabeloffset...
                    *wcos(azdir*(w.aztick(i)+w.azoffset));
            case 'out'
                pos(1)=pos(1)+w.azticklabeloffset...
                    *wsin(azdir*(w.aztick(i)+w.azoffset));
                pos(2)=pos(2)+w.azticklabeloffset...
                    *wcos(azdir*(w.aztick(i)+w.azoffset));
        end
        set(azh(i),'position',pos);
        set(azh(i),'units','data');
        % fix matlab bug
        pos=get(azh(i),'position');
        pos(3)=w.z;
        set(azh(i),'position',pos);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_minorticks(ax,w)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_minorgrid(ax,w)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wpv,lpv]=wedge_splitpv(w,pv)
% property names
wfields=fieldnames(w);

% extract linespec
if(mod(numel(pv),2))
    linespec=pv(1);
    pv(1)=[];
else
    linespec={};
end

% demand properties are strings
if(~iscellstr(pv(1:2:end)))
    error('seizmo:wedge:badInput',...
        'Properties must be strings!');
end

% separate
pv(1:2:end)=lower(pv(1:2:end));
widx=find(ismember(pv(1:2:end),wfields));
wpv=pv([2*widx-1; 2*widx]);
pv([2*widx-1 2*widx])=[];
lpv=[linespec pv];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wpv]=wedge_check(wpv)
% check cmdline wedge properties
for i=1:2:numel(wpv)
    switch lower(wpv{i})
        case {'axeslinewidth' 'fontsize' 'gridlinewidth' ...
              'minorgridlinewidth' 'ticklength'}
            if(~isreal(wpv{i+1}) || ~isscalar(wpv{i+1}) || wpv{i+1}<0)
                error('seizmo:wedge:badInput',...
                    '%s must be a positive real scalar!',upper(wpv{i}));
            end
        case {'axis' 'azgrid' 'azminorgrid' 'azminortick' ...
              'aztickvisible' 'azticklabelvisible' 'box' 'clipping' ...
              'plotbox' 'rgrid' 'rminorgrid' 'rminortick' ...
              'rtickvisible' 'rticklabelvisible'}
            checkstring(wpv{i:i+1},{'on' 'off'});
        case {'azaxislocation' 'tickdir'}
            checkstring(wpv{i:i+1},{'in' 'out'});
        case {'azcolor' 'backgroundcolor' 'gridcolor' 'minorgridcolor' ...
              'rcolor'}
            if(isnumeric(wpv{i+1}))
                if(~isequal(size(wpv{i+1}),[1 3]))
                    error('seizmo:wedge:badInput',...
                        'Numeric %s must be a RGB triplet!',upper(wpv{i}));
                end
            elseif(ischar(wpv{i+1}))
                if(~any(ismember(wpv{i+1},'ymcrgbwk')) ...
                        && ~strcmpi(wpv{i+1},'none'))
                    error('seizmo:wedge:badInput',...
                        'String %s must be a valid colorspec!',...
                        upper(wpv{i}));
                end
            else
                error('seizmo:wedge:badInput',...
                    '%s not valid!',upper(wpv{i}));
            end
        case {'azdir' 'raxislocation'}
            checkstring(wpv{i:i+1},{'cw' 'ccw'});
        case 'azlim'
            if(~isreal(wpv{i+1}) || ~isequal(size(wpv{i+1}),[1 2]) ...
                    || wpv{i+1}(1)>=wpv{i+1}(2))
                error('seizmo:wedge:badInput',...
                    'AZLIM must be a 1x2 array of [AZMIN AZMAX]!');
            end
        case {'azoffset' 'aztickdelta' 'azticklabeloffset' ...
              'rtickdelta' 'rticklabelcircaz' 'rticklabeloffset' 'z'}
            if(~isreal(wpv{i+1}) || ~isscalar(wpv{i+1}))
                error('seizmo:wedge:badInput',...
                    '%s must be a real-valued scalar!',upper(wpv{i}));
            end
        case {'aztick' 'rtick'}
            if(~isreal(wpv{i+1}) || ~isvector(wpv{i+1}) ...
                    || any(diff(wpv{i+1})<0))
                error('seizmo:wedge:badInput',...
                    ['%s must be a vector of monotonically ' ...
                    'increasing real values'],upper(wpv{i}));
            end
            wpv{i+1}=wpv{i+1}(:)'; % force row vector
        case {'azticklabel' 'rticklabel'}
            if(isempty(wpv{i+1}))
                wpv{i+1}=[];
            elseif(iscellstr(wpv{i+1}))
                % wonderful
            elseif(ischar(wpv{i+1}) && ndims(wpv{i+1})<=2)
                if(isvector(wpv{i+1}))
                    brk=[0 find(wpv{i+1}(:)'==124) numel(wpv{i+1})+1];
                    tmp=cell(numel(brk)-1,1);
                    for j=1:numel(brk)-1
                        tmp{j}=wpv{i+1}(brk(j)+1:brk(j+1)-1);
                    end
                    wpv{i+1}=tmp;
                else
                    wpv{i+1}=cellstr(wpv{i+1});
                end
            elseif(isnumeric(wpv{i+1}))
                wpv{i+1}=cellstr(num2str(wpv{i+1}(:)));
            else
                error('seizmo:wedge:badInput',...
                    '%s input not understood!',upper(wpv{i}));
            end
        case {'azticklabelhalign' 'rticklabelhalign'}
            checkstring(wpv{i:i+1},{'center' 'left' 'right' 'auto'});
        case {'azticklabelvalign' 'rticklabelvalign'}
            checkstring(wpv{i:i+1},...
                {'middle' 'top' 'cap' 'baseline' 'bottom'});
        case 'azunits'
            checkstring(wpv{i:i+1},{'degrees' 'radians'});
        case 'fontangle'
            checkstring(wpv{i:i+1},{'normal' 'italic' 'oblique'});
        case {'fontname' 'rticklabelunits' 'tag'}
            if(~ischar(wpv{i+1}))
                error('seizmo:wedge:badInput',...
                    '%s must be a string!',upper(wpv{i}));
            end
        case 'fontunits'
            checkstring(wpv{i:i+1},...
                {'points' 'normalized' 'inches' 'centimeters' 'pixels'});
        case 'fontweight'
            checkstring(wpv{i:i+1},{'normal' 'bold' 'light' 'demi'});
        case {'gridlinestyle' 'minorgridlinestyle' }
            checkstring(wpv{i:i+1},{'-' '--' ':' '-.' 'none'});
        case 'npts'
            if(~isreal(wpv{i+1}) || ~isscalar(wpv{i+1}) || wpv{i+1}<2)
                error('seizmo:wedge:badInput',...
                    'NPTS must be a positive real-valued scalar >1!');
            end
        case 'placement'
            checkstring(wpv{i:i+1},{'center' 'origin'});
        case 'rlim'
            if(~isreal(wpv{i+1}) || ~isequal(size(wpv{i+1}),[1 2]) ...
                    || wpv{i+1}(1)>=wpv{i+1}(2) || any(wpv{i+1}<0))
                error('seizmo:wedge:badInput',...
                    'RLIM must be a 1x2 positive array of [RMIN RMAX]!');
            end
        otherwise
            % Something I missed or changed?
            error('seizmo:wedge:badInput',...
                'Unknown option: %s !',wpv{i});
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=checkstring(name,str,valid)
valid=upper(valid);
if(~ischar(str) || ~any(strcmpi(str,valid)))
    error('seizmo:wedge:badInput',...
        ['%s must be one of the following strings:\n' ...
        sprintf('%s ',valid{:})],upper(name));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s]=appendpv2struct(s,pv)
for i=1:2:numel(pv); s.(pv{i})=pv{i+1}; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_update(ax,wpv)
error('sorry update not operative yet');
for i=1:2:numel(wpv)
    if(isequal(wpv{i+1},getappdata(ax,wpv{i}))); continue; end
    switch lower(wpv{i})
        case 'box'
        case 'clipping'
        case 'backgroundcolor'
        case 'fontangle'
        case 'fontname'
        case 'fontsize'
        case 'fontunits'
        case 'fontweight'
        case 'gridlinestyle'
        case 'axeslinewidth'
        case 'minorgridlinestyle'
        case 'tickdir'
        case 'ticklength'
        case 'visible'
        case 'raxislocation'
        case 'azaxislocation'
        case 'rcolor'
        case 'azcolor'
        case 'rlim'
        case 'rminorgrid'
        case 'rminortick'
        case 'rtick'
        case 'rticklabel'
        case 'azlim'
        case 'azminorgrid'
        case 'azminortick'
        case 'aztick'
        case 'azticklabel'
        case 'azdir'
        case 'azoffset'
        case 'azunits'
        case 'placement'
        case 'tag'
        case 'z'
        case 'npts'
        otherwise
            % Something I missed or changed?
            error('seizmo:wedge:badInput',...
                'Unknown option: %s !',wpv{i});
    end
    setappdata(ax,wpv{i},wpv{i+1});
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pv]=struct2pv(s)
pv={}; fields=fieldnames(s);
for i=1:numel(fields); pv=[pv {fields{i} s.(fields{i})}]; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,delta]=magicticks(limits,axwidth,fontsize)
%limits   - data limits
%axwidth  - axis width (normalized units)
%fontsize - fontsize (normalized units)
maxdigits=3.5;
valid=[1 2 2.5 5];
nmax=floor(axwidth/(fontsize*(maxdigits)))+1;
delta0=diff(limits);
delta=delta0;
cnt=1;
while(true)
    vi=mod(4-cnt,4)+1;
    exp=floor(log10(delta0))-floor((cnt-1)/4);
    delta=valid(vi)*10^exp;
    ticks=delta*(ceil(limits(1)/delta):floor(limits(2)/delta));
    if(any(numel(ticks)<[2 nmax])); cnt=cnt+1;
    else break;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,delta]=magicticks360(limits,axwidth,fontsize)
%limits   - data limits
%axwidth  - axis width (normalized units)
%fontsize - fontsize (normalized units)
maxdigits=3.5;
valid360=[1 2 2.5 5 10 15 20 30 40 45 60 90 120];
nmax=floor(axwidth/(fontsize*(maxdigits)))+1;
delta0=diff(limits);
delta=delta0;
for vi=13:-1:1
    delta=valid360(vi);
    ticks=delta*(ceil(limits(1)/delta):floor(limits(2)/delta));
    if(~any(numel(ticks)<[2 nmax])); return; end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=require_wedge(ax)
if(isempty(ax) || ~isscalar(ax) || ~ishandle(ax) ...
        || ~strcmp(get(ax,'createfcn'),'wedge'))
    error('seizmo:wedge:notWedge',...
        'No wedge plot found!');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%
function [s]=onoffstate(varargin)
% childonoff=onoffstate(parentonoff,childonoff)
for i=1:nargin
    if(strcmpi(varargin{i},'off')); s=varargin{i}; return; end
    s=varargin{end};
end
end

