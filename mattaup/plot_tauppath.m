function [varargout]=plot_tauppath(tt,varargin)
%PLOT_TAUPPATH    Plots tauppath output
%
%    Usage:    plot_tauppath(tt)
%              plot_tauppath(tt,'property',value,...)
%              ax=plot_tauppath(...)
%
%    Description:
%     PLOT_TAUPPATH(TT) plots the raypaths in TT on a new plot.  TT is the
%     output from TAUPPATH.  See that function for more details on
%     generating raypaths.
%
%     PLOT_TAUPPATH(TT,'PROPERTY',VALUE,...) allows some up-front plot
%     manipulation by adjusting specific properties.  Any unknown
%     properties are passed on to PLOT.  Valid properties:
%      EVANGLE    - angle to place event at (default is -45)
%      RAYCOLOR   - coloring of rays (colormap string, rgb triplet, name)
%                   (default is 'hsv')
%      STSYMBOL   - station marker symbol (default is '^')
%      STSIZE     - station marker size (default is 8)
%      STCOLOR    - station marker facecolor (default is 'r')
%      EVSYMBOL   - event marker symbol (default is 'p')
%      EVSIZE     - event marker size (default is 10)
%      EVCOLOR    - event marker facecolor (default is 'y')
%      BGCOLOR    - background color (rgb triplet, name) (default is 'k')
%      FGCOLOR    - foreground color (rgb triplet, name) (default is 'w')
%      SHOWLEGEND - true or false (default is true)
%      PARENT     - handle of axes to draw in (1 or 2 handles) (none)
%      DRAWEARTH  - draw earth when axes given to 'PARENT' option (false)
%
%     AX=PLOT_TAUPPATH(...) returns the handle to the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % Draw a variety of phases from an event to 3 random stations:
%     ax=plot_tauppath(tauppath,'raycolor','r','legend',false);
%     plot_tauppath(tauppath,'raycolor','g','parent',ax,'legend',false);
%     plot_tauppath(tauppath,'raycolor','b','parent',ax,'legend',false);
%
%    See also: TAUPPATH, PLOT_TAUPCURVE

%     Version History:
%        May  21, 2011 - initial version
%        Dec.  6, 2011 - fix raycolor bug
%        Feb. 24, 2012 - added lots of new properties
%        May   3, 2012 - tagged the grid, minor doc fix
%        Jan. 17, 2014 - added drawearth option, use parent for axes option
%        Jan. 27, 2014 - fix legend calls for octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 17:15 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>1 && ~mod(nargin,2))
    error('TauP:plot_tauppath:badNumOptions',...
        'Unpaired option(s)!');
end

% check tauppath output
if(~isstruct(tt) || any(~isfield(tt,...
        {'modelname' 'depth' 'distance' 'mindistance' 'phase' ...
        'puristphase' 'time' 'rayparameter' 'path'})))
    error('TauP:plot_tauppath:badInput',...
        'TT was not created by TAUPPATH!');
elseif(isempty(tt))
    error('TauP:plot_tauppath:badInput',...
        'TT is empty!');
end

% number of phases
nph=numel(tt);

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('TauP:plot_tauppath:badInput',...
        'One or more properties appear to not be given as strings!');
end

% strip axes arguments
[ax,varargin]=axparse(varargin{:});

% append defaults
varargin=[{'fgc' [] 'bgc' [] 'legend' true ...
    'evsym' 'p' 'evsize' 10 'evc' 'y' ...
    'stsym' '^' 'stsize' 8 'stc' 'r' ...
    'evang' -45 'rayc' 'hsv' 'de' false} varargin];

% extract pertinent parameters
delete=false(numel(varargin),1);
for i=1:2:numel(varargin)
    p=varargin{i};
    v=varargin{i+1};
    switch p
        case {'drawearth' 'drawe' 'de' 'earth' 'draw'}
            if(~isscalar(v) || ~islogical(v))
                error('TauP:plot_tauppath:badInput',...
                    'DRAWEARTH must be TRUE or FALSE!');
            end
            opt.drawearth=v;
            delete(i:i+1)=true;
        case {'evsymbol' 'evsym'}
            if(~ischar(v) || ~isscalar(v) ...
                    || ~any(lower(v)=='.ox+*sdv^<>ph'))
                error('TauP:plot_tauppath:badInput',...
                    'EVSYMBOL must be a valid single character symbol!');
            end
            opt.evsymbol=v;
            delete(i:i+1)=true;
        case {'evsize' 'evsz'}
            if(~isscalar(v) || ~isreal(v) || v<0)
                error('TauP:plot_tauppath:badInput',...
                    'EVSIZE must be a positive real-valued scalar!');
            end
            opt.evsize=v;
            delete(i:i+1)=true;
        case {'evcolor' 'evc'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_tauppath:badInput',...
                    'EVCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_tauppath:badInput',...
                    'EVCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.evcolor=v;
            delete(i:i+1)=true;
        case {'stsymbol' 'stsym'}
            if(~ischar(v) || ~isscalar(v) ...
                    || ~any(lower(v)=='.ox+*sdv^<>ph'))
                error('TauP:plot_tauppath:badInput',...
                    'STSYMBOL must be a valid single character symbol!');
            end
            opt.stsymbol=v;
            delete(i:i+1)=true;
        case {'stsize' 'stsz'}
            if(~isscalar(v) || ~isreal(v) || v<0)
                error('TauP:plot_tauppath:badInput',...
                    'STSIZE must be a positive real-valued scalar!');
            end
            opt.stsize=v;
            delete(i:i+1)=true;
        case {'stcolor' 'stc'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_tauppath:badInput',...
                    'STCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_tauppath:badInput',...
                    'STCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.stcolor=v;
            delete(i:i+1)=true;
        case {'fg' 'fgc' 'fgcolor'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_tauppath:badInput',...
                    'FGCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_tauppath:badInput',...
                    'FGCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.fgc=v;
            delete(i:i+1)=true;
        case {'bg' 'bgc' 'bgcolor'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_tauppath:badInput',...
                    'BGCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_tauppath:badInput',...
                    'BGCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.bgc=v;
            delete(i:i+1)=true;
        case {'showlegend' 'legend' 'leg'}
            if(~isscalar(v) || ~islogical(v))
                error('TauP:plot_tauppath:badInput',...
                    'SHOWLEGEND must be TRUE or FALSE!');
            end
            opt.legend=v;
            delete(i:i+1)=true;
        case {'evang' 'evangle'}
            if(~isscalar(v) || ~isreal(v))
                error('TauP:plot_tauppath:badInput',...
                    'EVANGLE must be a real-valued scalar!');
            end
            opt.evangle=v;
            delete(i:i+1)=true;
        case {'rayc' 'raycolor'}
            % try to decipher input
            if(ischar(v))
                try
                    tmp=str2func(v);
                    tmp=tmp(nph);
                    if(isequal(size(tmp),[nph 3]))
                        opt.raycolor=tmp;
                    else
                        error('TauP:plot_tauppath:nothing',...
                            'Ignore this!');
                    end
                catch
                    try
                        tmp=name2rgb(v);
                        opt.raycolor=tmp(ones(nph,1),:);
                    catch
                        error('TauP:plot_tauppath:badInput',...
                            'Could not decipher RAYCOLOR input!');
                    end
                end
            else
                if(isequal(size(v),[1 3]))
                    opt.raycolor=v(ones(nph,1),:);
                elseif(isequal(size(v),[nph 3]))
                    opt.raycolor=v;
                else
                    error('TauP:plot_tauppath:badInput',...
                        'Could not decipher RAYCOLOR input!');
                end
            end
            delete(i:i+1)=true;
    end
end
varargin(delete)=[];

% fix fg/bg
if(ischar(opt.bgc)); opt.bgc=name2rgb(opt.bgc); end
if(ischar(opt.fgc)); opt.fgc=name2rgb(opt.fgc); end
if(isempty(opt.fgc))
    if(isempty(opt.bgc))
        opt.fgc='w'; opt.bgc='k';
    else
        opt.fgc=invertcolor(opt.bgc,true);
    end
elseif(isempty(opt.bgc))
    opt.bgc=invertcolor(opt.fgc,true);
end

% new plot?
if(isempty(ax) || opt.drawearth)
    if(isempty(ax))
        % new figure if no parent
        fh=figure('color',opt.bgc,'defaulttextcolor',opt.fgc,...
            'defaultaxesxcolor',opt.fgc,'defaultaxesycolor',opt.fgc,...
            'name','TauP Ray Paths');
        ax=axes('parent',fh);
        new=true;
    else
        new=false;
    end
    
    % get hold state, set hold on
    isheld=ishold(ax);
    hold(ax,'on');
    
    % plot grid
    [cx,cy]=circle(6871);
    plot(ax,cx,cy,'color',[0.2 0.2 0.2],'linewidth',2,...
        'tag','plot_tauppath_grid');
    if(new); set(ax,'position',[0.025 0.05 0.95 0.9]); end
    [cx1,cy1]=circle(6871,180);
    [cx2,cy2]=circle(6771,180);
    plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',1,...
        'tag','plot_tauppath_grid');
    [cx1,cy1]=circle(6871,36);
    [cx2,cy2]=circle(6671,36);
    plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',2,...
        'tag','plot_tauppath_grid');
    [cx,cy]=circle(6371,4);
    plot(ax,cx(1:2:3),cy(1:2:3),'color',[0.2 0.2 0.2],'linewidth',2,...
        'tag','plot_tauppath_grid');
    plot(ax,cx(2:2:4),cy(2:2:4),'color',[0.2 0.2 0.2],'linewidth',2,...
        'tag','plot_tauppath_grid');
    
    % plot major discontinuities
    for i=[6371 3480 1220]
        [cx,cy]=circle(i);
        plot(ax,cx,cy,'color',opt.fgc,'linewidth',2,...
            'tag','plot_tauppath_major_discon');
    end
    
    % plot minor discontinuities
    for i=[5961 5711 3780]
        [cx,cy]=circle(i);
        plot(ax,cx,cy,'color',[0.5 0.5 0.5],'linewidth',1,...
            'tag','plot_tauppath_minor_discon');
    end
    axis(ax,'off','equal');
else
    % get hold state, set hold on
    isheld=ishold(ax);
    hold(ax,'on');
end

% loop over phases
for i=1:nph
    % plot phase path
    cx=(6371-tt(i).path.depth)...
        .*sin((opt.evangle+tt(i).path.distance)/180*pi);
    cy=(6371-tt(i).path.depth)...
        .*cos((opt.evangle+tt(i).path.distance)/180*pi);
    plot(ax,cx,cy,'color',opt.raycolor(i,:),...
        'displayname',tt(i).phase,'tag','tauppath_raypath',varargin{:});
end

% set default earthquake,station location
stloc=unique([tt.distance]);
stloc=stloc(stloc<=180);
cx=6371*sin((opt.evangle+[0 stloc])/180*pi);
cy=6371*cos((opt.evangle+[0 stloc])/180*pi);

% plot event and station markers
if(ischar(opt.evcolor)); opt.evcolor=name2rgb(opt.evcolor); end
if(ischar(opt.stcolor)); opt.stcolor=name2rgb(opt.stcolor); end
h1=plot(ax,cx(1),cy(1),opt.evsymbol,'markeredgecolor',opt.evcolor,...
    'markersize',opt.evsize,'markerfacecolor',opt.evcolor,...
    'tag','event_marker','displayname','event');
for i=2:numel(cx)
    h2=plot(ax,cx(i),cy(i),opt.stsymbol,'markeredgecolor',opt.stcolor,...
        'markersize',opt.stsize,'markerfacecolor',opt.stcolor,...
        'tag','station_marker','displayname','station');
end

% turn off hold if was off before
if(~isheld); hold(ax,'off'); end

% update legend
if(opt.legend)
    lh=legend(ax,[h1 h2 flipud(findobj(ax,'tag','tauppath_raypath'))'],...
        get([h1 h2 flipud(findobj(ax,'tag','tauppath_raypath'))'],...
        'displayname'),'location','westoutside');
    set(lh,'color','none','fontsize',6,'interpreter','none',...
        'edgecolor',opt.fgc,'textcolor',opt.fgc);
end

% output handles if desired
if(nargout); varargout{1}=ax; end

end
