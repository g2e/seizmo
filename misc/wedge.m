function [varargout]=wedge(varargin)
%WEDGE    Polar coordinate wedge plot
%
%    Usage:    wedge(azimuth,radius)
%              wedge(azimuth,radius,linespec)
%              wedge(azimuth,radius,image)
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
%     WEDGE(AZIMUTH,RADIUS,IMAGE) passes inputs to pcolor as a way to plot
%     images in polar coordinates.  AZIMUTH & RADIUS must be vectors or
%     matrices of size NRxNAZ where NR is the number of radii & NAZ is the
%     number of azimuths.  Sampling in azimuth and radius must be regular.
%     IMAGE must be NRxNAZ in size.  Inputs are automatically adjusted so
%     that pcolor plots the image at the correct positions.  Note that the
%     pixels do not curve and so the image has inherent inaccuracy in pixel
%     position that is exacerbated when pixels span a large azimuth range.
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
%      azminorgrid        'off'     % azimuthal minor grid
%      azminortick        'off'     % azimuthal minor ticks
%      azoffset           0         % 0 is vert, affected by azdir, azunits
%      aztick             []        % suitable values found by default
%      aztickdelta        []        % suitable delta found by default
%      azticklabel        []        % uses aztick values
%      azticklabelhalign  'auto'    % center/left/right/auto
%      azticklabeloffset  10        % radial offset in pixels
%      azticklabelvalign  'auto'    % middle/top/cap/baseline/bottom/auto
%      azticklabelvisible 'on'      % visibility of azimuthal tick labels
%      aztickmode         'auto'    % 'auto' or 'manual'
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
%      minorgridcolor     'k'       % color of minor grid lines
%      minorgridlinestyle ':'       % linestyle of minor grid lines
%      minorgridlinewidth 0.5       % linewidth of minor grid lines
%      npts               51        % number of points in azimuthal axis
%      placement          'center'  % center/origin
%      plotbox            'off'     % on/off (clip wedge during pan/zoom)
%      raxislocation      'cw'      % cw/ccw (only drawn for wedge)
%      rcolor             'k'       % color of radial axis+ticks+labels
%      rgrid              'on'      % visibility of radial grid lines
%      rlim               []        % defaults to [0 max(radius)]
%      rminorgrid         'off'     % azimuthal minor grid
%      rminortick         'off'     % azimuthal minor ticks
%      rtick              []        % suitable values found by default
%      rtickdelta         []        % suitable values found by default
%      rticklabel         []        % uses rtick values
%      rticklabelcircaz   10        % r axis label azimuth when circle
%      rticklabelhalign   'auto'    % center/left/right/auto
%      rticklabeloffset   10        % radial tick label offset in pixels
%      rticklabelvalign   'auto'    % middle/top/cap/baseline/bottom/auto
%      rticklabelvisible  'on'      % visibility of radial tick labels
%      rticklabelunits    ''        % units string to append to labels
%      rtickmode          'auto'    % 'auto' or 'manual'
%      rtickvisible       'on'      % visibility of radial ticks
%      tickdir            'in'      % direction of ticks
%      ticklength         0.02      % fraction of rlim difference
%      wedgetag           ''        % tag for wedge axes
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
%     [VAL1,VAL2,...]=WEDGE(AX,{'PROP1' 'PROP2' ...}) returns multiple
%     properties of the wedge axes.
%
%    Notes:
%     - The handle of the axes is the parent object of the plotted object.
%       Thus: AX=get(H,'parent');
%
%    Examples:
%     % Polar plot:
%     wedge(360*rand(10,1),rand(10,1),'azdir','ccw','azoffset',-90);
%
%     % peaks plot:
%     [x,y,z]=peaks;
%     wedge(x*120,y,z,'azlim',[-180 180]);
%
%     % Radioactive wedge:
%     wedge(axes,[],[],'azlim',[-15 15],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%     wedge(axes,[],[],'azlim',[105 135],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%     wedge(axes,[],[],'azlim',[225 255],'rlim',[0.5 1.5],...
%           'placement','origin','backgroundcolor','r')
%
%    See also: WGRID, POLAR, ROSE, COMPASS

%     Version History:
%        Feb. 22, 2011 - initial version
%        Feb. 23, 2011 - fix missing handle output, add padding for text,
%                        changed halign & offset defaults
%        Feb. 24, 2011 - pcolor option
%        Apr. 27, 2012 - fix minor breakage in p/v checking
%        May   1, 2012 - updating many parameters is now available
%        May   2, 2012 - many more enhancements
%        Sep. 12, 2012 - doc update: expect image, not pcolor input
%        Jan. 27, 2014 - use axparse instead of axescheck for octave,
%                        adjust newplot call only for octave
%        Mar. 17, 2014 - clipping (hidden opt & image only) defaults to off
%        Mar. 18, 2014 - image with <300 pixels now behind axes/grid
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2014 at 18:35 GMT

% todo
% - drawing
%   - improve clipping (need to figure this one out)
%   - updating tick length/dir should update label position
%   - labels should be a constant distance from ticks (callback)
%   - low res image pixels don't curve making the plot inaccurate
% - updating wedges
%   - anything but a complete redraw right now
%     - redraw requires some work
%       - how to respect the plotted objects when doing this
% - rlabel/azlabel (separate functions)
%   - where to position when circle
%   - this will need callbacks for positioning and more
% - multiple wedges in the same axes
%   - for cool stuff like princeton hotspot plots
% - wedge3 (separate code)
%   - cylindrical vs spherical system
%     - wedge3 should be spherical while wedge is cylindrical
%       - wedge rotation would draw axes extension in the cylindrical
%         domain...
% - deg:min:sec stuff
% - something akin to plotyy (degrees vs radians, radius vs depth)

% find axes input
[ax,varargin]=axparse(varargin{:});
nargs=numel(varargin);

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
                && ~isequal(size(varargin{1}),size(varargin{2}),...
                size(varargin{3})))
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
            [ax,w]=wedge_create(ax,varargin{2},w);
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
        wedge_update(ax,wedge_check(varargin));
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
wedge_cleanup(ax,held);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    degstep=x(1,2)-x(1,1);
    radstep=y(2,1)-y(1,1);
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

% <300 pcolor block fix
if(prod(size(x)-1)<300)
    movekids(h,'back');
    movekids(findall(ax,'tag','wedge_bgpatch'),'back');
end

% cleanup
wedge_cleanup(ax,held);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax,w]=wedge_create(ax,y,w)
%WEDGE_CREATE    Creates wedge plot

% clear plot & set hold to on
if(ishandle(ax)); hold(ax,'off'); end
if(exist('OCTAVE_VERSION','builtin')==5)
    newplot;
else
    ax=newplot(ax);
end
hold(ax,'on');

% tag the underlying axes
set(ax,'tag',w.wedgetag);

% set font
set(ax,'fontname',w.fontname,'fontunits',w.fontunits,...
    'fontsize',w.fontsize,'fontweight',w.fontweight,...
    'fontangle',w.fontangle);

% first determine auto rlim & tlim
w=wedge_auto_rlim(w,y);
w=wedge_auto_azlim(w);

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
w.azminorgrid='off';
w.azminortick='off';
w.azoffset=0; % 0 is vertical, affected by azdir, azunits
w.aztick=[];
w.aztickdelta=[];
w.azticklabel=[];
w.azticklabelhalign='auto'; % center/left/right/auto
w.azticklabeloffset=10; % outward offset in pixels
w.azticklabelvalign='auto'; % middle/top/cap/baseline/bottom/auto
w.azticklabelvisible='on';
w.aztickmode='auto'; % auto/manual
w.aztickvisible='on';
w.azunits='degrees'; % degrees/radians
w.backgroundcolor='w';
w.box='on';
w.clipping='off'; % only works (sorta) for image
w.fontangle=get(0,'defaulttextfontangle');
w.fontname=get(0,'defaulttextfontname');
w.fontsize=get(0,'defaulttextfontsize');
w.fontunits=get(0,'defaulttextfontunits');
w.fontweight=get(0,'defaulttextfontweight');
w.gridcolor='k';
w.gridlinestyle=':';
w.gridlinewidth=0.5;
w.minorgridcolor='k';
w.minorgridlinestyle=':';
w.minorgridlinewidth=0.5;
w.npts=51;
w.placement='center'; % center/origin
w.plotbox='off'; % on/off (clip wedge during pan/zoom)
w.raxislocation='cw'; % cw/ccw for wedge
w.rcolor='k';
w.rgrid='on';
w.rlim=[];
w.rminorgrid='off';
w.rminortick='off';
w.rtick=[];
w.rtickdelta=[];
w.rticklabel=[];
w.rticklabelcircaz=10; % r axis label azimuth when circle
w.rticklabelhalign='auto'; % center/left/right/auto
w.rticklabeloffset=10; % outward offset in pixels
w.rticklabelvalign='auto'; % middle/top/cap/baseline/bottom/auto
w.rticklabelvisible='on';
w.rticklabelunits='';
w.rtickmode='auto'; % auto/manual
w.rtickvisible='on';
w.tickdir='in';
w.ticklength=0.02; % normalized to rlim difference
w.wedgetag='';
w.z=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_axes(ax,w)

% get metaparameters
m=wedge_metaparameters(w);

% axes lines
axesinx=[];  axesiny=[];  % not necessarily drawn
axescwx=[];  axescwy=[];  % not necessarily drawn
axesccwx=[]; axesccwy=[]; % not necessarily drawn
axesoutx=m.wsin(m.azdir*...
    (linspace(w.azlim(1),w.azlim(2),w.npts)'+w.azoffset));
axesouty=m.wcos(m.azdir*...
    (linspace(w.azlim(1),w.azlim(2),w.npts)'+w.azoffset));
if(m.hollow)
    % hollow, so draw the hole
    axesinx=w.rlim(1)*axesoutx; % note that positions above were
    axesiny=w.rlim(1)*axesouty; % normalized to make this easy
end
axesoutx=w.rlim(2)*axesoutx; % now place in
axesouty=w.rlim(2)*axesouty; % correct position
if(~m.circ)
    % not a circle so we need axes extending in the radial direction
    axescwx=[w.rlim(1)*m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset));...
             w.rlim(2)*m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    axescwy=[w.rlim(1)*m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset));...
             w.rlim(2)*m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    axesccwx=[w.rlim(1)*m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset));...
              w.rlim(2)*m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
    axesccwy=[w.rlim(1)*m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset));...
              w.rlim(2)*m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
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
    line(axesinx,axesiny,w.z(ones(numel(axesinx),1)),...
        'parent',ax,'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_in',...
        'clipping',w.plotbox,'visible',axisstate('in',w));
end
line(axesoutx,axesouty,w.z(ones(numel(axesoutx),1)),...
    'parent',ax,'color',w.azcolor,'linewidth',w.axeslinewidth,...
    'handlevisibility','off','tag','wedge_axis_out',...
    'clipping',w.plotbox,'visible',axisstate('out',w));
if(~isempty(axescwx))
    line(axescwx,axescwy,w.z(ones(numel(axescwx),1)),...
        'parent',ax,'color',w.rcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_cw',...
        'clipping',w.plotbox,'visible',axisstate('cw',w));
end
if(~isempty(axesccwx))
    line(axesccwx,axesccwy,w.z(ones(numel(axesccwx),1)),...
        'parent',ax,'color',w.rcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_axis_ccw',...
        'clipping',w.plotbox,'visible',axisstate('ccw',w));
end

% axes placement + 15% padding
axis(ax,'equal');
if(strcmp(w.placement,'origin'))
    axis(ax,w.rlim(2)*[-1.15 1.15 -1.15 1.15]);
else % center
    axis(ax,'tight'); % shifts wedge to center
    xrng=get(ax,'xlim');
    yrng=get(ax,'ylim');
    set(ax,'xlim',xrng+.15*[-1 1]*diff(xrng)/2,...
           'ylim',yrng+.15*[-1 1]*diff(yrng)/2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w]=wedge_ticks(ax,w)

% get metaparameters
m=wedge_metaparameters(w);

% auto tick setup
[rwidth,azwidth,fsize]=wedge_widths(ax,w,m);

% get auto ticks
switch w.aztickmode
    % ticks => mode=manual, delta ignored
    % delta => mode=auto, ticks replaced
    % auto => ticks replaced, delta kept
    % manual => ticks kept, delta ignored
    case 'auto'
        % ignore ticks given, use delta if given
        if(isempty(w.aztickdelta))
            % full auto
            if(m.azdeg) % degrees
                if(m.circ)
                    w.aztick=magicticks360(w.azlim,azwidth,fsize);
                else
                    w.aztick=magicticks(w.azlim,azwidth,fsize);
                end
            else % radians
                w.aztick=magicticks(w.azlim./pi,azwidth,fsize);
                w.aztick=w.aztick*pi;
            end
            % don't have 2 ticks at the same spot
            if(m.circ)
                if(m.azdeg)
                    if(w.aztick(1)+360==w.aztick(end))
                        w.aztick(end)=[];
                    end
                else % radians
                    if(w.aztick(1)+2*pi==w.aztick(end));
                        w.aztick(end)=[];
                    end
                end
            end
        else % delta given
            w.aztick=w.aztickdelta*(ceil(...
                w.azlim(1)/w.aztickdelta):floor(w.azlim(2)/w.aztickdelta));
            if(m.circ)
                if(m.azdeg)
                    if(w.aztick(1)+360==w.aztick(end))
                        w.aztick(end)=[];
                    end
                else % radians
                    if(w.aztick(1)+2*pi==w.aztick(end))
                        w.aztick(end)=[];
                    end
                end
            end
        end
    case 'manual'
        % use ticks given
end
switch w.rtickmode
    case 'auto'
        % ignore ticks given, use delta if given
        if(isempty(w.rtickdelta))
            w.rtick=magicticks(w.rlim,rwidth,fsize);
            if(~w.rtick(1)); w.rtick(1)=[]; end
        else
            w.rtick=w.rtickdelta*(ceil(...
                w.rlim(1)/w.rtickdelta):floor(w.rlim(2)/w.rtickdelta));
            if(w.rtick(1)==0); w.rtick(1)=[]; end
        end
    case 'manual'
        % use ticks given
end

% tick mark length
abstlen=diff(w.rlim)*w.ticklength;

% tick lines
tickinx=[];  tickiny=[];
tickoutx=[]; tickouty=[];
tickcwx=[];  tickcwy=[];
tickccwx=[]; tickccwy=[];
if(~isempty(w.aztick))
    tickoutr=w.rlim(2)+[zeros(1,numel(w.aztick));
        -m.tdir*abstlen(ones(1,numel(w.aztick)))];
    tickoutx=tickoutr.*m.wsin(m.azdir*(w.aztick([1 1],:)+w.azoffset));
    tickouty=tickoutr.*m.wcos(m.azdir*(w.aztick([1 1],:)+w.azoffset));
    if(m.hollow)
        % inner axes exists so ticks too
        tickinr=w.rlim(1)+[zeros(1,numel(w.aztick));
            m.tdir*abstlen(ones(1,numel(w.aztick)))];
        tickinx=tickinr.*m.wsin(m.azdir*(w.aztick([1 1],:)+w.azoffset));
        tickiny=tickinr.*m.wcos(m.azdir*(w.aztick([1 1],:)+w.azoffset));
    end
end
if(~isempty(w.rtick) && ~m.circ)
    % cw/ccw axes exist so ticks too
    tickcwx=w.rtick*m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset));
    tickcwy=w.rtick*m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset));
    tickcwx=[tickcwx; tickcwx+m.tdir*abstlen...
        *m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    tickcwy=[tickcwy; tickcwy-m.tdir*abstlen...
        *m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    tickccwx=w.rtick*m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset));
    tickccwy=w.rtick*m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset));
    tickccwx=[tickccwx; tickccwx-m.tdir*abstlen...
        *m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
    tickccwy=[tickccwy; tickccwy+m.tdir*abstlen...
        *m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
end

% draw ticks
if(~isempty(tickoutx))
    line(tickoutx,tickouty,w.z(ones(2,numel(w.aztick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_out',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('out',w),w.aztickvisible));
end
if(~isempty(tickinx))
    line(tickinx,tickiny,w.z(ones(2,numel(w.aztick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_in',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('in',w),w.aztickvisible));
end
if(~isempty(tickcwx))
    line(tickcwx,tickcwy,w.z(ones(2,numel(w.rtick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_cw',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('cw',w),w.rtickvisible));
end
if(~isempty(tickccwx))
    line(tickccwx,tickccwy,w.z(ones(2,numel(w.rtick))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_tick_ccw',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('ccw',w),w.rtickvisible));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_grid(ax,w)

% get metaparameters
m=wedge_metaparameters(w);

% skip if no radial ticks
if(~isempty(w.rtick))
    % radial grid (circles)
    %x=r.*wsin(azdir*(az+w.azoffset))
    %y=r.*wcos(azdir*(az+w.azoffset))
    azpnts=linspace(w.azlim(1),w.azlim(2),w.npts)';
    rgridx=w.rtick(ones(w.npts,1),:) ...
        .*m.wsin(m.azdir*(azpnts(:,ones(1,numel(w.rtick)))+w.azoffset));
    rgridy=w.rtick(ones(w.npts,1),:) ...
        .*m.wcos(m.azdir*(azpnts(:,ones(1,numel(w.rtick)))+w.azoffset));
    
    % plot radial grid
    line(rgridx,rgridy,w.z(ones(w.npts,numel(w.rtick))),'parent',ax,...
        'color',w.gridcolor,'linewidth',w.gridlinewidth,...
        'linestyle',w.gridlinestyle,'handlevisibility','off',...
        'visible',onoffstate(w.axis,w.rgrid),...
        'tag','wedge_grid_r','clipping',w.plotbox);
end

% skip if no azimuthal ticks
if(~isempty(w.aztick))
    % azimuthal grid (spokes)
    rpnts=w.rlim';
    azgridx=rpnts(:,ones(1,numel(w.aztick))) ...
        .*m.wsin(m.azdir*(w.aztick(ones(2,1),:)+w.azoffset));
    azgridy=rpnts(:,ones(1,numel(w.aztick))) ...
        .*m.wcos(m.azdir*(w.aztick(ones(2,1),:)+w.azoffset));
    
    % plot azimuthal grid
    line(azgridx,azgridy,w.z(ones(2,numel(w.aztick))),'parent',ax,...
        'color',w.gridcolor,'linewidth',w.gridlinewidth,...
        'linestyle',w.gridlinestyle,'handlevisibility','off',...
        'visible',onoffstate(w.axis,w.azgrid),...
        'tag','wedge_grid_az','clipping',w.plotbox);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_ticklabels(ax,w)

% get metaparameters
m=wedge_metaparameters(w);

% auto labels
if(isempty(w.azticklabel))
    if(m.azdeg)
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
if(m.circ)
    % labels are at position specified by rticklabelcircaz
    rlabelx=w.rtick.*m.wsin(m.azdir*(w.rticklabelcircaz+w.azoffset));
    rlabely=w.rtick.*m.wcos(m.azdir*(w.rticklabelcircaz+w.azoffset));
else
    switch w.raxislocation
        case 'cw'
            rlabelx=w.rtick.*m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset));
            rlabely=w.rtick.*m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset));
        case 'ccw'
            rlabelx=w.rtick.*m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset));
            rlabely=w.rtick.*m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset));
    end
end

% azimuthal tick label positioning
switch w.azaxislocation
    case 'in'
        if(m.hollow)
            azlabelx=w.rlim(1)*m.wsin(m.azdir*(w.aztick+w.azoffset));
            azlabely=w.rlim(1)*m.wcos(m.azdir*(w.aztick+w.azoffset));
        else % do not plot labels
            azlabelx=[];
            azlabely=[];
        end
    case 'out'
        azlabelx=w.rlim(2)*m.wsin(m.azdir*(w.aztick+w.azoffset));
        azlabely=w.rlim(2)*m.wcos(m.azdir*(w.aztick+w.azoffset));
end

% auto horizontal alignment
rhalign=cellstr(w.rticklabelhalign);
rhalign=rhalign(ones(numel(w.rtick),1));
if(strcmp(w.rticklabelhalign,'auto'))
    if(m.circ)
        rhalign(:)={'center'};
    else
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

% auto vertical alignment
rvalign=cellstr(w.rticklabelvalign);
rvalign=rvalign(ones(numel(w.rtick),1));
if(strcmp(w.rticklabelvalign,'auto'))
    rvalign(:)={'middle'}; % radial autos to middle for now
end
azvalign=cellstr(w.azticklabelvalign);
azvalign=azvalign(ones(numel(w.aztick),1));
if(strcmp(w.azticklabelvalign,'auto') && ~isempty(azlabelx))
    azvalign(:)={'middle'}; % azimuthal autos to middle for now
    %{
    % auto pushes text outwards from plot
    % out & +y => bottom
    %        0 => middle
    %       -y => top
    % in &  +y => top
    %        0 => middle
    %       -y => bottom
    for i=1:numel(azlabelx)
        switch w.azaxislocation
            case 'in'
                if(azlabely(i)>eps('single'))
                    azvalign{i}='top';
                elseif(azlabely(i)<-eps('single'))
                    azvalign{i}='bottom';
                else
                    azvalign{i}='middle';
                end
            case 'out'
                if(azlabely(i)>eps('single'))
                    azvalign{i}='bottom';
                elseif(azlabely(i)<-eps('single'))
                    azvalign{i}='top';
                else
                    azvalign{i}='middle';
                end
        end
    end
    %}
end

% plot labels
rh=nan(numel(w.rtick),1);
for i=1:numel(w.rtick)
    w.rticklabel{i}=[w.rticklabel{i} w.rticklabelunits];
    rh(i)=text(rlabelx(i),rlabely(i),w.z,w.rticklabel{i},'parent',ax,...
        'horizontalalignment',rhalign{i},'clipping',w.plotbox,...
        'visible',onoffstate(w.axis,w.rticklabelvisible),...
        'verticalalignment',rvalign{i},'color',w.rcolor,...
        'handlevisibility','off','tag','wedge_rlabels',...
        'fontunits',w.fontunits,'fontsize',w.fontsize,...
        'fontname',w.fontname,'fontangle',w.fontangle,...
        'fontweight',w.fontweight);
end
if(~isempty(azlabelx))
    azh=nan(numel(w.aztick),1);
    for i=1:numel(w.aztick)
        if(~m.azdeg && strcmp(w.azticklabel{i},'0'))
            % no units added on
        elseif(~m.azdeg && strcmp(w.azticklabel{i},'1'))
            % drop the 1
            w.azticklabel{i}=m.azticklabelunits;
        else
            w.azticklabel{i}=[w.azticklabel{i} m.azticklabelunits];
        end
        azh(i)=text(azlabelx(i),azlabely(i),w.z,w.azticklabel{i},...
            'parent',ax,'horizontalalignment',azhalign{i},...
            'clipping',w.plotbox,'handlevisibility','off',...
            'visible',onoffstate(w.axis,w.azticklabelvisible),...
            'verticalalignment',azvalign{i},'color',w.azcolor,...
            'tag','wedge_azlabels','fontunits',w.fontunits,...
            'fontsize',w.fontsize,'fontname',w.fontname,...
            'fontangle',w.fontangle,'fontweight',w.fontweight);
    end
end

% offset labels
for i=1:numel(w.rtick)
    set(rh(i),'units','pixels');
    pos=get(rh(i),'position');
    % modify position
    if(m.circ) % radially if circular
        pos(1)=pos(1)+w.rticklabeloffset...
            *m.wsin(m.azdir*(w.rticklabelcircaz+w.azoffset));
        pos(2)=pos(2)+w.rticklabeloffset...
            *m.wcos(m.azdir*(w.rticklabelcircaz+w.azoffset));
    else % azimuthally otherwise
        switch w.raxislocation
            case 'cw'
                pos(1)=pos(1)-w.rticklabeloffset...
                    *m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset));
                pos(2)=pos(2)+w.rticklabeloffset...
                    *m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset));
            case 'ccw'
                pos(1)=pos(1)+w.rticklabeloffset...
                    *m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset));
                pos(2)=pos(2)-w.rticklabeloffset...
                    *m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset));
        end
    end
    set(rh(i),'position',pos);
    set(rh(i),'units','data');
    % fix matlab bug
    pos=get(rh(i),'position');
    pos(3)=w.z;
    set(rh(i),'position',pos);
end
if(~isempty(azlabelx))
    for i=1:numel(w.aztick)
        set(azh(i),'units','pixels');
        pos=get(azh(i),'position');
        % modify position
        switch w.azaxislocation
            case 'in'
                pos(1)=pos(1)-w.azticklabeloffset...
                    *m.wsin(m.azdir*(w.aztick(i)+w.azoffset));
                pos(2)=pos(2)-w.azticklabeloffset...
                    *m.wcos(m.azdir*(w.aztick(i)+w.azoffset));
            case 'out'
                pos(1)=pos(1)+w.azticklabeloffset...
                    *m.wsin(m.azdir*(w.aztick(i)+w.azoffset));
                pos(2)=pos(2)+w.azticklabeloffset...
                    *m.wcos(m.azdir*(w.aztick(i)+w.azoffset));
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

% get metaparameters
m=wedge_metaparameters(w);

% auto tick setup
[rwidth,azwidth,fsize]=wedge_widths(ax,w,m);

% get auto minorticks
azminorticks=[];
rminorticks=[];
if(~isempty(w.aztick))
    if(m.azdeg)
        azminorticks=magicminorticks(w.azlim,azwidth,fsize,w.aztick);
    else % radians
        azminorticks=magicminorticks(w.azlim./pi,azwidth,fsize,w.aztick);
        azminorticks=azminorticks.*pi;
    end
end
if(~isempty(w.rtick))
    if(isscalar(w.rtick))
        rminorticks=magicminorticks(w.rlim,rwidth,fsize,[0 w.rtick]);
    else
        rminorticks=magicminorticks(w.rlim,rwidth,fsize,w.rtick);
    end
end

% minortick mark length is 1/2 normal tick mark length
abstlen=diff(w.rlim)*w.ticklength/2;

% tick lines
tickinx=[];  tickiny=[];
tickoutx=[]; tickouty=[];
tickcwx=[];  tickcwy=[];
tickccwx=[]; tickccwy=[];
if(~isempty(azminorticks))
    tickoutr=w.rlim(2)+[zeros(1,numel(azminorticks));
        -m.tdir*abstlen(ones(1,numel(azminorticks)))];
    tickoutx=tickoutr.*m.wsin(m.azdir*(azminorticks([1 1],:)+w.azoffset));
    tickouty=tickoutr.*m.wcos(m.azdir*(azminorticks([1 1],:)+w.azoffset));
    if(m.hollow)
        % inner axes exists so ticks too
        tickinr=w.rlim(1)+[zeros(1,numel(azminorticks));
            m.tdir*abstlen(ones(1,numel(azminorticks)))];
        tickinx=tickinr.*m.wsin(...
            m.azdir*(azminorticks([1 1],:)+w.azoffset));
        tickiny=tickinr.*m.wcos(...
            m.azdir*(azminorticks([1 1],:)+w.azoffset));
    end
end
if(~isempty(rminorticks) && ~m.circ)
    % cw/ccw axes exist so ticks too
    tickcwx=rminorticks*m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset));
    tickcwy=rminorticks*m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset));
    tickcwx=[tickcwx; tickcwx+m.tdir*abstlen...
        *m.wcos(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    tickcwy=[tickcwy; tickcwy-m.tdir*abstlen...
        *m.wsin(m.azdir*(w.azlim(m.cw)+w.azoffset))];
    tickccwx=rminorticks*m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset));
    tickccwy=rminorticks*m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset));
    tickccwx=[tickccwx; tickccwx-m.tdir*abstlen...
        *m.wcos(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
    tickccwy=[tickccwy; tickccwy+m.tdir*abstlen...
        *m.wsin(m.azdir*(w.azlim(m.ccw)+w.azoffset))];
end

% draw ticks
if(~isempty(tickoutx))
    line(tickoutx,tickouty,w.z(ones(2,numel(azminorticks))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_minortick_out',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('out',w),w.azminortick));
end
if(~isempty(tickinx))
    line(tickinx,tickiny,w.z(ones(2,numel(azminorticks))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_minortick_in',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('in',w),w.azminortick));
end
if(~isempty(tickcwx))
    line(tickcwx,tickcwy,w.z(ones(2,numel(rminorticks))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_minortick_cw',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('cw',w),w.rminortick));
end
if(~isempty(tickccwx))
    line(tickccwx,tickccwy,w.z(ones(2,numel(rminorticks))),'parent',ax,...
        'color',w.azcolor,'linewidth',w.axeslinewidth,...
        'handlevisibility','off','tag','wedge_minortick_ccw',...
        'clipping',w.plotbox,...
        'visible',onoffstate(axisstate('ccw',w),w.rminortick));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_minorgrid(ax,w)

% get metaparameters
m=wedge_metaparameters(w);

% auto tick setup
[rwidth,azwidth,fsize]=wedge_widths(ax,w,m);

% get auto minorticks
azminorticks=[];
rminorticks=[];
if(~isempty(w.aztick))
    if(m.azdeg)
        azminorticks=magicminorticks(w.azlim,azwidth,fsize,w.aztick);
    else % radians
        azminorticks=magicminorticks(w.azlim./pi,azwidth,fsize,w.aztick);
        azminorticks=azminorticks.*pi;
    end
end
if(~isempty(w.rtick))
    if(isscalar(w.rtick))
        rminorticks=magicminorticks(w.rlim,rwidth,fsize,[0 w.rtick]);
    else
        rminorticks=magicminorticks(w.rlim,rwidth,fsize,w.rtick);
    end
end

% skip if no radial minorticks
if(~isempty(rminorticks))
    % radial grid (circles)
    %x=r.*wsin(azdir*(az+w.azoffset))
    %y=r.*wcos(azdir*(az+w.azoffset))
    azpnts=linspace(w.azlim(1),w.azlim(2),w.npts)';
    rgridx=rminorticks(ones(w.npts,1),:).*m.wsin(...
        m.azdir*(azpnts(:,ones(1,numel(rminorticks)))+w.azoffset));
    rgridy=rminorticks(ones(w.npts,1),:).*m.wcos(...
        m.azdir*(azpnts(:,ones(1,numel(rminorticks)))+w.azoffset));
    
    % plot radial grid
    line(rgridx,rgridy,w.z(ones(w.npts,numel(rminorticks))),'parent',ax,...
        'color',w.minorgridcolor,'linewidth',w.minorgridlinewidth,...
        'linestyle',w.minorgridlinestyle,'handlevisibility','off',...
        'visible',onoffstate(w.axis,w.rminorgrid),...
        'tag','wedge_minorgrid_r','clipping',w.plotbox);
end

% skip if no azimuthal ticks
if(~isempty(azminorticks))
    % azimuthal grid (spokes)
    rpnts=w.rlim';
    azgridx=rpnts(:,ones(1,numel(azminorticks))) ...
        .*m.wsin(m.azdir*(azminorticks(ones(2,1),:)+w.azoffset));
    azgridy=rpnts(:,ones(1,numel(azminorticks))) ...
        .*m.wcos(m.azdir*(azminorticks(ones(2,1),:)+w.azoffset));
    
    % plot azimuthal grid
    line(azgridx,azgridy,w.z(ones(2,numel(azminorticks))),'parent',ax,...
        'color',w.minorgridcolor,'linewidth',w.minorgridlinewidth,...
        'linestyle',w.minorgridlinestyle,'handlevisibility','off',...
        'visible',onoffstate(w.axis,w.azminorgrid),...
        'tag','wedge_minorgrid_az','clipping',w.plotbox);
end

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
i=1;
while(i<numel(wpv))
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
        case {'aztickmode' 'rtickmode'}
            checkstring(wpv{i:i+1},{'auto' 'manual'});
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
        case {'azoffset' 'azticklabeloffset' ...
              'rticklabelcircaz' 'rticklabeloffset' 'z'}
            if(~isreal(wpv{i+1}) || ~isscalar(wpv{i+1}))
                error('seizmo:wedge:badInput',...
                    '%s must be a real-valued scalar!',upper(wpv{i}));
            end
        case {'aztickdelta' 'rtickdelta'}
            if(~isreal(wpv{i+1}) || ~isscalar(wpv{i+1}))
                error('seizmo:wedge:badInput',...
                    '%s must be a real-valued scalar!',upper(wpv{i}));
            end
            % requires ticks drawn 2 times
            wpv=[wpv(1:i-1) {[wpv{i}(1:end-4) 'tickmode'] 'auto'} ...
                wpv(i:end)];
            i=i+2;
        case {'aztick' 'rtick'}
            if(~isreal(wpv{i+1}) || ~isvector(wpv{i+1}) ...
                    || any(diff(wpv{i+1})<0))
                error('seizmo:wedge:badInput',...
                    ['%s must be a vector of monotonically ' ...
                    'increasing real values'],upper(wpv{i}));
            end
            wpv{i+1}=wpv{i+1}(:)'; % force row vector
            % requires ticks drawn 2 times
            wpv=[wpv(1:i-1) {[wpv{i}(1:end-4) 'tickmode'] 'manual'} ...
                wpv(i:end)];
            i=i+2;
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
                {'middle' 'top' 'cap' 'baseline' 'bottom' 'auto'});
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
    i=i+2;
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

% struct of current parameters
w=getappdata(ax);

% loop over each p/v pair
for i=1:2:numel(wpv)
    % skip any parameters that are unchanged
    if(isequal(wpv{i+1},w.(wpv{i}))); continue; end
    
    % assign to struct
    w.(wpv{i})=wpv{i+1};
    
    % operation depends on parameter
    switch lower(wpv{i})
        case 'box'
            % change box visibility
            wedge_redo_box(ax,w);
        case 'backgroundcolor'
            set(findall(ax,'tag','wedge_bgpatch'),'facecolor',wpv{i+1});
        case 'aztickvisible'
            % tick visibility (respect box too)
            set(findall(ax,'tag','wedge_tick_in'),...
                'visible',onoffstate(axisstate('in',w),wpv{i+1}));
            set(findall(ax,'tag','wedge_tick_out'),...
                'visible',onoffstate(axisstate('out',w),wpv{i+1}));
        case 'rtickvisible'
            % tick visibility (respect box too)
            set(findall(ax,'tag','wedge_tick_cw'),...
                'visible',onoffstate(axisstate('cw',w),wpv{i+1}));
            set(findall(ax,'tag','wedge_tick_ccw'),...
                'visible',onoffstate(axisstate('ccw',w),wpv{i+1}));
        case {'rticklabelvisible' 'azticklabelvisible'}
            % tick label visibility
            set(findall(ax,'tag',['wedge_' wpv{i}(1:end-16) 'labels']),...
                'visible',onoffstate(w.axis,wpv{i+1}));
        case {'rgrid' 'azgrid'}
            % grid visibility
            set(findall(ax,'tag',['wedge_grid_' wpv{i}(1:end-4)]),...
                'visible',onoffstate(w.axis,wpv{i+1}));
        case {'gridlinestyle' 'gridlinewidth' 'gridcolor'}
            % changes grid only
            delete(findall(ax,'-regexp','tag','^wedge_grid_'));
            wedge_grid(ax,w);
        case {'rminorgrid' 'azminorgrid'}
            % minor grid visibility
            set(findall(ax,'tag',['wedge_minorgrid_' wpv{i}(1:end-4)]),...
                'visible',onoffstate(w.axis,wpv{i+1}));
        case {'minorgridlinestyle' 'minorgridlinewidth' 'minorgridcolor'}
            % changes minor grid only
            delete(findall(ax,'-regexp','tag','^wedge_minorgrid_'));
            wedge_minorgrid(ax,w);
        case 'axeslinewidth'
            % line width of axes & ticks
            set(findall(ax,'-regexp','tag','^wedge_axis'),...
                'linewidth',wpv{i+1});
            set(findall(ax,'-regexp','tag','^wedge_tick'),...
                'linewidth',wpv{i+1});
        case {'tickdir' 'ticklength'}
            % redraw ticks
            delete(findall(ax,'-regexp','tag','^wedge_tick_'));
            wedge_ticks(ax,w);
        case 'axis'
            % hide/show axes, ticks, labels, background, grid
            set(findall(ax,'tag','wedge_axis_in'),...
                'visible',axisstate('in',w));
            set(findall(ax,'tag','wedge_axis_out'),...
                'visible',axisstate('out',w));
            set(findall(ax,'tag','wedge_axis_cw'),...
                'visible',axisstate('cw',w));
            set(findall(ax,'tag','wedge_axis_ccw'),...
                'visible',axisstate('ccw',w));
            set(findall(ax,'tag','wedge_tick_cw'),...
                'visible',onoffstate(axisstate('cw',w),w.rtickvisible));
            set(findall(ax,'tag','wedge_tick_ccw'),...
                'visible',onoffstate(axisstate('ccw',w),w.rtickvisible));
            set(findall(ax,'tag','wedge_tick_in'),...
                'visible',onoffstate(axisstate('in',w),w.aztickvisible));
            set(findall(ax,'tag','wedge_tick_out'),...
                'visible',onoffstate(axisstate('out',w),w.aztickvisible));
            set(findall(ax,'tag','wedge_rlabels'),...
                'visible',onoffstate(wpv{i+1},w.rticklabelvisible));
            set(findall(ax,'tag','wedge_azlabels'),...
                'visible',onoffstate(wpv{i+1},w.azticklabelvisible));
            set(findall(ax,'tag','wedge_bgpatch'),'visible',wpv{i+1});
            set(findall(ax,'tag','wedge_grid_r'),...
                'visible',onoffstate(wpv{i+1},w.rgrid));
            set(findall(ax,'tag','wedge_grid_az'),...
                'visible',onoffstate(wpv{i+1},w.azgrid));
            set(findall(ax,'tag','wedge_minorgrid_r'),...
                'visible',onoffstate(wpv{i+1},w.rminorgrid));
            set(findall(ax,'tag','wedge_minorgrid_az'),...
                'visible',onoffstate(wpv{i+1},w.azminorgrid));
            set(findall(ax,'tag','wedge_minortick_cw'),'visible',...
                onoffstate(axisstate('cw',w),w.rminortick));
            set(findall(ax,'tag','wedge_minortick_ccw'),'visible',...
                onoffstate(axisstate('ccw',w),w.rminortick));
            set(findall(ax,'tag','wedge_minortick_in'),'visible',...
                onoffstate(axisstate('in',w),w.azminortick));
            set(findall(ax,'tag','wedge_minortick_out'),'visible',...
                onoffstate(axisstate('out',w),w.azminortick));
        case {'raxislocation' 'azaxislocation'}
            % redraw labels
            delete(findall(ax,'-regexp','tag','^wedge_(r|az)labels$'));
            wedge_ticklabels(ax,w);
            % change box visibility
            wedge_redo_box(ax,w);
        case 'rcolor'
            % color of axes, ticks, labels
            set(findall(ax,'tag','wedge_axis_cw'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_axis_ccw'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_tick_cw'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_tick_ccw'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_rlabels'),'color',wpv{i+1});
        case 'azcolor'
            % color of axes, ticks, labels
            set(findall(ax,'tag','wedge_axis_in'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_axis_out'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_tick_in'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_tick_out'),'color',wpv{i+1});
            set(findall(ax,'tag','wedge_azlabels'),'color',wpv{i+1});
        case 'rminortick'
            % minortick visibility (respect box too)
            set(findall(ax,'tag','wedge_minortick_cw'),...
                'visible',onoffstate(axisstate('cw',w),wpv{i+1}));
            set(findall(ax,'tag','wedge_minortick_ccw'),...
                'visible',onoffstate(axisstate('ccw',w),wpv{i+1}));
        case 'azminortick'
            % minortick visibility (respect box too)
            set(findall(ax,'tag','wedge_minortick_in'),...
                'visible',onoffstate(axisstate('in',w),wpv{i+1}));
            set(findall(ax,'tag','wedge_minortick_out'),...
                'visible',onoffstate(axisstate('out',w),wpv{i+1}));
        case {'rtick' 'aztick' 'rtickdelta' 'aztickdelta' 'rtickmode' ...
                'aztickmode' 'fontangle' 'fontname' 'fontsize' ...
                'fontunits' 'fontweight'}
            % redo ticks fully
            wedge_redo_ticks_full(ax,w,wpv,i);
        case {'azticklabel' 'azticklabelhalign' 'azticklabeloffset' ...
                'azticklabelvalign' 'rticklabel' 'rticklabelhalign' ...
                'rticklabeloffset' 'rticklabelvalign' ...
                'rticklabelcircaz' 'rticklabelunits'}
            % redraw labels
            delete(findall(ax,'-regexp','tag','^wedge_(r|az)labels$'));
            wedge_ticklabels(ax,w);
        case 'placement'
            % axes placement + 15% padding
            axis(ax,'equal');
            if(strcmp(w.placement,'origin'))
                axis(ax,w.rlim(2)*[-1.15 1.15 -1.15 1.15]);
            else % center
                % this has trouble for data not clipped
                axis(ax,'tight'); % shifts wedge to center
                xrng=get(ax,'xlim');
                yrng=get(ax,'ylim');
                set(ax,'xlim',xrng+.15*[-1 1]*diff(xrng)/2,...
                       'ylim',yrng+.15*[-1 1]*diff(yrng)/2);
            end
            % redo ticks fully
            wedge_redo_ticks_full(ax,w,wpv,i);
        case 'wedgetag'
            % tag the underlying axes
            set(ax,'tag',wpv{i+1});
        case {'npts' 'z' 'azunits' 'azoffset' 'azdir' ...
                'azlim' 'rlim' 'clipping'}
            % Something unimplemented
            error('seizmo:wedge:unimplementedInput',...
                'Update not allowed for option: %s !',wpv{i});
        otherwise
            % Something I missed or changed?
            error('seizmo:wedge:badInput',...
                'Unknown option: %s !',wpv{i});
    end
    setappdata(ax,wpv{i},wpv{i+1});
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_redo_ticks_full(ax,w,wpv,i)

% redraw ticks, tick labels, grid
delete(findall(ax,'-regexp','tag','^wedge_tick_'));
delete(findall(ax,'-regexp','tag','^wedge_(r|az)labels$'));
delete(findall(ax,'-regexp','tag','^wedge_grid_'));
delete(findall(ax,'-regexp','tag','^wedge_minortick_'));
delete(findall(ax,'-regexp','tag','^wedge_minorgrid_'));
if(strcmpi(wpv{i},'rtickdelta')); w.rtick=[]; end
if(strcmpi(wpv{i},'aztickdelta')); w.aztick=[]; end
w=wedge_ticks(ax,w);
wedge_ticklabels(ax,w);
wedge_grid(ax,w);
wedge_minorticks(ax,w);
wedge_minorgrid(ax,w);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_redo_box(ax,w)

% hide/show axes & ticks
set(findall(ax,'tag','wedge_axis_in'),...
    'visible',axisstate('in',w));
set(findall(ax,'tag','wedge_axis_out'),...
    'visible',axisstate('out',w));
set(findall(ax,'tag','wedge_axis_cw'),...
    'visible',axisstate('cw',w));
set(findall(ax,'tag','wedge_axis_ccw'),...
    'visible',axisstate('ccw',w));
set(findall(ax,'tag','wedge_tick_cw'),...
    'visible',onoffstate(axisstate('cw',w),w.rtickvisible));
set(findall(ax,'tag','wedge_tick_ccw'),...
    'visible',onoffstate(axisstate('ccw',w),w.rtickvisible));
set(findall(ax,'tag','wedge_tick_in'),...
    'visible',onoffstate(axisstate('in',w),w.aztickvisible));
set(findall(ax,'tag','wedge_tick_out'),...
    'visible',onoffstate(axisstate('out',w),w.aztickvisible));
set(findall(ax,'tag','wedge_minortick_cw'),'visible',...
    onoffstate(axisstate('cw',w),w.rminortick));
set(findall(ax,'tag','wedge_minortick_ccw'),'visible',...
    onoffstate(axisstate('ccw',w),w.rminortick));
set(findall(ax,'tag','wedge_minortick_in'),'visible',...
    onoffstate(axisstate('in',w),w.azminortick));
set(findall(ax,'tag','wedge_minortick_out'),'visible',...
    onoffstate(axisstate('out',w),w.azminortick));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pv]=struct2pv(s)
pv={}; fields=fieldnames(s);
for i=1:numel(fields); pv=[pv {fields{i} s.(fields{i})}]; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w]=wedge_auto_rlim(w,y)
if(isempty(w.rlim))
    rmax=max(abs(y(~isinf(y))));
    if(isempty(rmax)); rmax=1; end
    w.rlim=[0 rmax];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w]=wedge_auto_azlim(w)
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
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m]=wedge_metaparameters(w)

% setup for azimuth units & directions
m.tdir=1; if(strcmp(w.tickdir,'out')); m.tdir=-1; end
m.azaxpo=w.rlim(1);
if(strcmp(w.azaxislocation,'out'))
    m.azaxpo=w.rlim(2);
end
m.azdir=1; if(strcmp(w.azdir,'ccw')); m.azdir=-1; end
m.cw=(3-m.azdir)/2; m.ccw=(3+m.azdir)/2; % cw/ccw = 1 or 2 (indices)
if(strcmp(w.azunits,'radians')); m.azdeg=false; m.wsin=@sin; m.wcos=@cos;
else m.azdeg=true; m.wsin=@sind; m.wcos=@cosd;
end

% label units & full circle
if(m.azdeg)
    m.azticklabelunits='^o';
    m.circ=diff(w.azlim)>=360;
else
    m.azticklabelunits='\pi';
    m.circ=diff(w.azlim)>=2*pi;
end

% hollow?
if(w.rlim(1)<=0); m.hollow=false; else m.hollow=true; end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rwidth,azwidth,fsize]=wedge_widths(ax,w,m)

% radial, azimuthal length
xlim=get(ax,'xlim');
units=get(ax,'units');
set(ax,'units','pixels');
%pos=get(ax,'position');
pos=true_axis_position(ax);
set(ax,'units',units);
rwidth=diff(w.rlim)*pos(3)/diff(xlim);
if(m.azdeg); azwidth=diff(w.azlim)*pi/180*m.azaxpo*pos(3)/diff(xlim);
else azwidth=diff(w.azlim)*m.azaxpo*pos(3)/diff(xlim);
end

% font size
set(ax,'fontunits','pixels');
fsize=get(ax,'fontsize');
set(ax,'fontunits',w.fontunits);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,delta]=magicticks(limits,axwidth,fontsize)
%limits   - data limits
%axwidth  - axis width (normalized units)
%fontsize - fontsize (normalized units)
maxdigits=4;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,delta]=magicticks360(limits,axwidth,fontsize)
%limits   - data limits
%axwidth  - axis width (normalized units)
%fontsize - fontsize (normalized units)
maxdigits=4;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,delta]=magicminorticks(limits,axwidth,fontsize,ticks1)
%limits   - data limits
%axwidth  - axis width (normalized units)
%fontsize - fontsize (normalized units)
maxdigits=4;
valid=[1 2 2.5 5];
delta0=diff(limits);
delta1=ticks1(2)-ticks1(1);
nmax=floor(axwidth*delta1/delta0/(fontsize*maxdigits))+1;
delta=delta1;
cnt=1;
while(true)
    vi=mod(4-cnt,4)+1;
    exp=floor(log10(delta1))-floor((cnt-1)/4);
    delta=valid(vi)*10^exp;
    tticks=delta*(ceil(ticks1(1)/delta):floor(ticks1(2)/delta));
    ticks=delta*(ceil(limits(1)/delta):floor(limits(2)/delta));
    if(any(numel(tticks)<[3 nmax])); cnt=cnt+1;
    else break;
    end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s]=onoffstate(varargin)
% childonoff=onoffstate(parentonoff,childonoff)
for i=1:nargin
    if(strcmpi(varargin{i},'off')); s=varargin{i}; return; end
    s=varargin{end};
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s]=axisstate(axis,w)
switch axis
    case 'in'
        switch w.azaxislocation
            case 'in'
                s=onoffstate(w.axis);
            case 'out'
                s=onoffstate(w.axis,w.box);
        end
    case 'out'
        switch w.azaxislocation
            case 'in'
                s=onoffstate(w.axis,w.box);
            case 'out'
                s=onoffstate(w.axis);
        end
    case 'cw'
        switch w.raxislocation
            case 'cw'
                s=onoffstate(w.axis);
            case 'ccw'
                s=onoffstate(w.axis,w.box);
        end
    case 'ccw'
        switch w.raxislocation
            case 'cw'
                s=onoffstate(w.axis,w.box);
            case 'ccw'
                s=onoffstate(w.axis);
        end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=wedge_cleanup(ax,held)

% "hidden" tag to identify wedge axes
set(ax,'createfcn','wedge');

% new plot or just adding to it?
if(~held)
    % restore hold state
    hold(ax,'off');
    
    % hide underlying axis but make the x/y labels visible
    axis(ax,'off');
    set(get(ax,'xlabel'),'visible','on');
    set(get(ax,'ylabel'),'visible','on');
end

end

