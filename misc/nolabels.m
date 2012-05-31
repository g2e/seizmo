function nolabels(ax,xy)
%NOLABELS    Removes tick and axis labels
%
%    Usage:    nolabels()
%              nolabels(ax)
%              nolabels(ax,option)
%
%    Description:
%     NOLABELS() removes x/y tick labels and axis labels for the current
%     axis.  This does not remove the ticks.
%
%     NOLABELS(AX) removes x/y tick labels and axis labels for the
%     specified axes.
%
%     NOLABELS(AX,OPTION) specifies the labels to be removed:
%      'xy'      -- (DEFAULT) removes x & y axes & tick labels
%      'x'       -- removes x axis & tick labels
%      'y'       -- removes y axis & tick labels
%      'xytick'  -- removes x & y axes tick labels
%      'xtick'   -- removes x axis tick labels
%      'ytick'   -- removes y axis tick labels
%      'xylabel' -- removes x & y axes labels
%      'xlabel'  -- removes x axis label
%      'ylabel'  -- removes y axis label
%
%    Notes:
%
%    Examples:
%     % make a figure with 4x3 arrangement of
%     % subplots and edit labels and ticks
%     figure;
%     ax=makesubplots(5,3,1:12);
%     ax=reshape(ax,3,4);
%     tax=ax';
%     axexpand(ax,15);
%     nolabels(tax(5:12),'y');
%     nolabels(ax(1:9),'x');
%     noticks(tax(5:12),'y');
%     noticks(ax(1:9),'x');
%
%    See also: NOTICKS, NOCOLORBARS, NOTITLES, SETFONTS, MAKESUBPLOTS,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              AXSTRETCH, AXEXPAND, AXMOVE

%     Version History:
%        Aug.  4, 2010 - initial version
%        Aug.  7, 2010 - added 'xtick' 'ytick' 'xytick' 'xlabel' 'ylabel' &
%                        'xylabel' options
%        Aug.  8, 2010 - doc update
%        May   9, 2012 - made options a bit more flexible
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   9, 2012 at 17:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% defaults
if(nargin<1 || isempty(ax)); ax=gca; end
if(nargin<2 || isempty(xy)); xy='xy'; end

% check inputs
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:nolabels:badInput',...
        'AX must be a valid axes handle!');
elseif(~ischar(xy))
    error('seizmo:nolabels:badInput',...
        'XY must be either ''x'' ''y'' or ''xy''!');
end

% implement on axis of choice
switch lower(xy)
    case 'xy'
        set(ax,'xticklabel',[],'yticklabel',[]);
        tmp=get(ax,'xlabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
        tmp=get(ax,'ylabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case 'x'
        set(ax,'xticklabel',[]);
        tmp=get(ax,'xlabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case 'y'
        set(ax,'yticklabel',[]);
        tmp=get(ax,'ylabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case {'xylabel' 'xylabels'}
        tmp=get(ax,'xlabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
        tmp=get(ax,'ylabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case {'xlabel' 'xlabels'}
        tmp=get(ax,'xlabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case {'ylabel' 'ylabels'}
        tmp=get(ax,'ylabel');
        if(iscell(tmp)); tmp=cat(1,tmp{:}); end
        set(tmp,'string',[]);
    case {'xytick' 'xyticks'}
        set(ax,'xticklabel',[],'yticklabel',[]);
    case {'xtick' 'xticks'}
        set(ax,'xticklabel',[]);
    case {'ytick' 'yticks'}
        set(ax,'yticklabel',[]);
    otherwise
        error('seizmo:nolabels:badInput',...
            'XY input unknown!');
end

end
