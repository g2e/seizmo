function [varargout]=plot_tauppath(tt,varargin)
%PLOT_TAUPPATH    Plots tauppath output
%
%    Usage:    plot_tauppath(tt)
%              plot_tauppath(tt,'option',value,...)
%              h=plot_tauppath(...)
%
%    Description:
%     PLOT_TAUPPATH(TT) plots the raypaths in TT on a new plot.  TT is the
%     output from TAUPPATH.  See that function for more details on
%     generating raypaths.
%
%     PLOT_TAUPPATH(TT,'OPTION',VALUE,...) passes option/value pairs to
%     PLOT.  This is useful for defining the axis ('parent') or the width
%     of the paths ('linewidth').  A couple extra parameters are available:
%      'evangle'  - angle to place event at (-45 is the default)
%      'raycolor' - coloring of rays (colormap string, rgb values, name)
%
%     H=PLOT_TAUPPATH(...) returns the handles to the raypaths.
%
%    Notes:
%
%    Examples:
%     % Draw a variety of phases from an event to 3 random stations:
%     plot_tauppath(tauppath,'raycolor','r')
%     plot_tauppath(tauppath,'raycolor','g','parent',gca)
%     plot_tauppath(tauppath,'raycolor','b','parent',gca)
%     legend off
%
%    See also: TAUPPATH, PLOT_TAUPCURVE

%     Version History:
%        May  21, 2011 - initial version
%        Dec.  6, 2011 - fix raycolor bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  6, 2011 at 17:15 GMT

% todo:
% - fgc/bgc
% - station marker options
% - event marker options

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>1 && ~mod(nargin,2))
    error('matTaup:plot_tauppath:badNumOptions',...
        'Unpaired option(s)!');
end

% check taupcurve output
if(~isstruct(tt) || any(~isfield(tt,...
        {'time' 'distance' 'depth' 'phase' 'rayparameter' 'path'})))
    error('matTaup:plot_tauppath:badInput',...
        'TT was not created by TAUPPATH!');
elseif(isempty(tt))
    error('matTaup:plot_tauppath:badInput',...
        'TT is empty!');
end

% number of phases
nph=numel(tt);

% was an axis handle given?
[ax,varargin]=axescheck(varargin{:});
if(isempty(ax)); newax=true; else newax=false; end

% append defaults
varargin=[{'evang' -45 'raycolor' 'hsv'} varargin];

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('matTaup:plot_tauppath:badInput',...
        'One or more parameters appear to not be given as strings!');
end

% extract pertinent parameters
% - event start location
% - raypath coloring
delete=false(numel(varargin),1);
for i=1:2:numel(varargin)
    p=varargin{i};
    v=varargin{i+1};
    switch p
        case {'evang' 'evangle'}
            P.evang=v;
            delete(i:i+1)=true;
        case {'rayc' 'raycolor'}
            % try to decipher input
            if(ischar(v))
                try
                    tmp=str2func(v);
                    tmp=tmp(nph);
                    if(isequal(size(tmp),[nph 3]))
                        P.raycolor=tmp;
                    else
                        error('matTaup:plot_tauppath:nothing',...
                            'Ignore this!');
                    end
                catch
                    try
                        tmp=name2rgb(v);
                        P.raycolor=tmp(ones(nph,1),:);
                    catch
                        error('matTaup:plot_tauppath:badInput',...
                            'Could not decipher RAYCOLOR input!');
                    end
                end
            else
                if(isequal(size(v),[1 3]))
                    P.raycolor=v(ones(nph,1),:);
                elseif(isequal(size(v),[nph 3]))
                    P.raycolor=v;
                else
                    error('matTaup:plot_tauppath:badInput',...
                        'Could not decipher RAYCOLOR input!');
                end
            end
            delete(i:i+1)=true;
    end
end
varargin(delete)=[];

% new plot?
if(newax)
    % new figure if no parent
    fh=figure('color','k','name','TauP Ray Paths');
    ax=axes('parent',fh);
    hold(ax,'on'); isheld=true;
    
    % plot grid
    [cx,cy]=circle(6871);
    plot(ax,cx,cy,'color',[0.2 0.2 0.2],'linewidth',2);
    set(ax,'position',[0.025 0.05 0.95 0.9]);
    [cx1,cy1]=circle(6871,180);
    [cx2,cy2]=circle(6771,180);
    plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',1);
    [cx1,cy1]=circle(6871,36);
    [cx2,cy2]=circle(6671,36);
    plot(ax,[cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',2);
    [cx,cy]=circle(6371,4);
    plot(ax,cx(1:2:3),cy(1:2:3),'color',[0.2 0.2 0.2],'linewidth',2);
    plot(ax,cx(2:2:4),cy(2:2:4),'color',[0.2 0.2 0.2],'linewidth',2);
    
    % plot major discontinuities
    for i=[6371 3480 1220]
        [cx,cy]=circle(i);
        plot(ax,cx,cy,'w','linewidth',2);
    end
    
    % plot minor discontinuities
    for i=[5961 5711 3780]
        [cx,cy]=circle(i);
        plot(ax,cx,cy,'color',[0.5 0.5 0.5],'linewidth',1);
    end
    axis(ax,'off','equal');
else
    % get hold state, set hold on
    isheld=ishold(ax);
    hold(ax,'on');
end

% loop over phases
ph=nan(1,nph); pn=cell(1,nph);
for i=1:nph
    % plot phase path
    cx=(6371-tt(i).path.depth).*sin((P.evang+tt(i).path.distance)/180*pi);
    cy=(6371-tt(i).path.depth).*cos((P.evang+tt(i).path.distance)/180*pi);
    ph(i)=plot(ax,cx,cy,'color',P.raycolor(i,:),...
        'displayname',tt(i).phase,'tag','tauppath_raypath',varargin{:});
    pn{i}=tt(i).phase;
end

% set default earthquake,station location
stloc=unique([tt.distance]);
stloc=stloc(stloc<=180);
cx=6371*sin((P.evang+[0 stloc])/180*pi);
cy=6371*cos((P.evang+[0 stloc])/180*pi);

% plot event and station markers
h1=plot(ax,cx(1),cy(1),'yp','markersize',10,'markerfacecolor','y',...
    'tag','event_marker','displayname','event');
for i=2:numel(cx)
    h2=plot(ax,cx(i),cy(i),'r^','markersize',8,'markerfacecolor','r',...
        'tag','station_marker','displayname','station');
end

% turn off hold if was off before
if(~isheld); hold(ax,'off'); end

% update legend
lh=legend([h1 h2 findobj(ax,'tag','tauppath_raypath')'],...
    'location','westoutside');
set(lh,'color','none','edgecolor','w','fontsize',6,'interpreter','none',...
    'textcolor','w');

% output handles if desired
if(nargout); varargout{1}=ph; end

end
