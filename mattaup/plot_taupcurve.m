function [varargout]=plot_taupcurve(tt,varargin)
%PLOT_TAUPCURVE    Plots taupcurve output
%
%    Usage:    plot_taupcurve(tt)
%              plot_taupcurve(tt,'property',value,...)
%           ax=plot_taupcurve(...)
%
%    Description:
%     PLOT_TAUPCURVE(TT) plots the output struct TT from TAUPCURVE in 2 new
%     figures.  The travel time curves vs distance are presented in one
%     figure and the rayparameter curves vs distance in another.
%
%     PLOT_TAUPCURVE(TT,'PROPERTY',VALUE,...) allows some up-front plot
%     manipulation by adjusting specific properties.  Any unknown
%     properties are passed on to PLOT.  Valid properties:
%      BGCOLOR    - background color (rgb triplet, name) (default is 'k')
%      FGCOLOR    - foreground color (rgb triplet, name) (default is 'w')
%      DRAW       - which plots to make ('both'=default, 'time', 'rayp')
%      SHOWLEGEND - true or false (default is true)
%      AXES       - handle(s) of axes to draw in (1 or 2 handles) (none)
%      CURVECOLOR - colormap, rgb triplet, name ('hsv')
%
%     AX=PLOT_TAUPCURVE(...) returns the handles to the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % A colorful example:
%     tt=taupcurve('z',500);
%     plot_taupcurve(tt,'bgcolor','purple','curvecolor','fire')
%
%    See also: PLOT_TAUPPATH, PLOT_TAUPCURVE_DT, TAUPCURVE

%     Version History:
%        Feb. 24, 2012 - initial version
%        Jan. 27, 2014 - fix legend calls for octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 17:15 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>1 && ~mod(nargin,2))
    error('TauP:plot_taupcurve:badNumOptions',...
        'Unpaired option(s)!');
end

% check taupcurve output
if(~isstruct(tt) || any(~isfield(tt,...
        {'modelname' 'depth' 'distance' 'mindistance' ...
        'phase' 'puristphase' 'time' 'rayparameter'})))
    error('TauP:plot_taupcurve:badInput',...
        'TT was not created by TAUPCURVE!');
elseif(isempty(tt))
    error('TauP:plot_taupcurve:badInput',...
        'TT is empty!');
end

% number of curves
nc=numel(tt);

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('TauP:plot_taupcurve:badInput',...
        'One or more properties appear to not be given as strings!');
end

% strip axes arguments
[ax,varargin]=axparse(varargin{:});

% append defaults
varargin=[{'fgc' [] 'bgc' [] 'draw' 'both' 'legend' true ...
    'curvecolor' 'hsv'} varargin];

% extract pertinent parameters
delete=false(numel(varargin),1);
for i=1:2:numel(varargin)
    p=varargin{i};
    v=varargin{i+1};
    switch p
        case {'fg' 'fgc' 'fgcolor'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_taupcurve:badInput',...
                    'FGCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_taupcurve:badInput',...
                    'FGCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.fgc=v;
            delete(i:i+1)=true;
        case {'bg' 'bgc' 'bgcolor'}
            if(~ischar(v) && ~isnumeric(v))
                error('TauP:plot_taupcurve:badInput',...
                    'BGCOLOR must be a valid colorstring or RGB triplet!');
            elseif(isnumeric(v) && ~isempty(v) && ~isequal(size(v),[1 3]))
                error('TauP:plot_taupcurve:badInput',...
                    'BGCOLOR must be a valid colorstring or RGB triplet!');
            end
            opt.bgc=v;
            delete(i:i+1)=true;
        case {'showlegend' 'legend' 'leg'}
            if(~isscalar(v) || ~islogical(v))
                error('TauP:plot_taupcurve:badInput',...
                    'SHOWLEGEND must be TRUE or FALSE!');
            end
            opt.legend=v;
            delete(i:i+1)=true;
        case {'draw'}
            if(isempty(v)); continue; end
            if(~ischar(v) || size(v,1)~=1 ...
                    || ~ismember(lower(v),{'time' 'rayp' 'both'}))
                error('TauP:plot_taupcurve:badInput',...
                    'DRAW must be ''TIME'' ''RAYP'' or ''BOTH''!');
            end
            opt.draw=v;
            delete(i:i+1)=true;
        case {'color' 'curvecolor' 'curve' 'ccolor' 'curvec' 'cc'}
            if(isempty(v)); continue; end
            % try to decipher input
            if(ischar(v))
                try
                    tmp=str2func(v);
                    tmp=tmp(nc);
                    if(isequal(size(tmp),[nc 3]))
                        opt.color=tmp;
                    else
                        error('TauP:plot_taupcurve:nothing',...
                            'Ignore this!');
                    end
                catch
                    try
                        tmp=name2rgb(v);
                        opt.color=tmp(ones(nc,1),:);
                    catch
                        error('TauP:plot_taupcurve:badInput',...
                            'Could not decipher CURVECOLOR input!');
                    end
                end
            else
                if(isequal(size(v),[1 3]))
                    opt.color=v(ones(nph,1),:);
                elseif(isequal(size(v),[nc 3]))
                    opt.color=v;
                else
                    error('TauP:plot_taupcurve:badInput',...
                        'Could not decipher CURVECOLOR input!');
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

% new plots or existing
if(strcmpi(opt.draw,'both'))
    if(numel(ax)==2); newplot=false; else newplot=true; end
else
    if(numel(ax)==1); newplot=false; else newplot=true; end
end

if(ismember(lower(opt.draw),{'time' 'both'}))
    % initialize dist vs time plot
    if(newplot)
        fh(1)=figure('color',opt.bgc,'defaulttextcolor',opt.fgc,...
            'defaultaxesxcolor',opt.fgc,'defaultaxesycolor',opt.fgc,...
            'name','TauP Travel Time Curves');
        pos=get(fh(1),'position');
        set(fh(1),'position',...
            [pos(1) pos(2)-pos(4)*.75 pos(3) pos(4)*1.75]);
        ax(1)=axes('parent',fh(1));
        set(ax(1),'color','none',...
            'yaxislocation','right','position',[0 0.1 0.9 0.8]);
    end
    title(ax(1),'TauP Travel Time Curves');
    xlabel(ax(1),'Distance (deg)');
    ylabel(ax(1),'Time (sec)');
    held=ishold(ax(1));
    if(~held); hold(ax(1),'on'); end
    box(ax(1),'on');
    
    % loop over phases
    for ii=1:nc
        plot(ax(1),tt(ii).mindistance,tt(ii).time,...
            'color',opt.color(ii,:),'displayname',tt(ii).phase,...
            'tag','taupcurve_timecurve',varargin{:});
    end
    
    % force distance range from 0 to 180
    xlim(ax(1),[0 180])
    
    % legend
    if(opt.legend)
        lh(1)=legend(ax(1),...
            flipud(findobj(ax(1),'tag','taupcurve_timecurve')),...
            get(flipud(findobj(ax(1),'tag','taupcurve_timecurve')),...
            'displayname'),'location','westoutside');
        set(lh(1),'color','none','edgecolor',opt.fgc,...
            'textcolor',opt.fgc,'fontsize',6,'interpreter','none');
    end
    if(~held); hold(ax(1),'off'); end
end

if(ismember(lower(opt.draw),{'rayp' 'both'}))
    % deal with handles
    if(strcmpi(opt.draw,'both')); np=2; else np=1; end
    
    % initialize dist vs ray parameter plot
    if(newplot)
        fh(np)=figure('color',opt.bgc,'defaulttextcolor',opt.fgc,...
            'defaultaxesxcolor',opt.fgc,'defaultaxesycolor',opt.fgc,...
            'name','TauP Ray Parameter Curves');
        pos=get(fh(np),'position');
        set(fh(np),'position',...
            [pos(1) pos(2)-pos(4)*.75 pos(3) pos(4)*1.75]);
        ax(np)=axes('parent',fh(np));
        set(ax(np),'color','none',...
            'yaxislocation','right','position',[0 0.1 0.9 0.8]);
    end
    title(ax(np),'TauP Ray Parameter Curves');
    xlabel(ax(np),'Distance (deg)');
    ylabel(ax(np),'Ray Parameter (sec/deg)');
    held=ishold(ax(np));
    if(~held); hold(ax(np),'on'); end
    box(ax(np),'on');
    
    % plot rayparameter vs dist
    for ii=1:nc
        plot(ax(np),tt(ii).mindistance,tt(ii).rayparameter,...
            'color',opt.color(ii,:),'displayname',tt(ii).phase,...
            'tag','taupcurve_rayparametercurve',varargin{:});
    end
    
    % force distance range from 0 to 180
    xlim(ax(np),[0 180])
    
    % legend
    if(opt.legend)
        lh(np)=legend(ax(np),...
            flipud(findobj(ax(np),'tag','taupcurve_rayparametercurve')),...
            get(flipud(findobj(ax(np),...
            'tag','taupcurve_rayparametercurve')),'displayname'),...
            'location','westoutside');
        set(lh(np),'color','none','edgecolor',opt.fgc,...
            'textcolor',opt.fgc,'fontsize',6,'interpreter','none');
    end
    if(~held); hold(ax(np),'off'); end
end

% output
if(nargout); varargout{1}=ax; end

end
