function nolabels(ax,xy)
%NOLABELS    Removes tick and axis labels
%
%    Usage:    nolabels()
%              nolabels(ax)
%              nolabels(ax,'x')
%              nolabels(ax,'y')
%
%    Description:
%     NOLABELS() removes x/y tick labels and axis labels for the current
%     axis.  This does not remove the ticks.
%
%     NOLABELS(AX) removes x/y tick labels and axis labels for the
%     specified axes.
%
%     NOLABELS(AX,'X') only removes the x-axis tick labels and x-axis
%     label.
%
%     NOLABELS(AX,'Y') only removes the y-axis tick labels and y-axis
%     label.
%
%    Notes:
%
%    Examples:
%     % make a figure with 4x3 arrangement of subplots, then expand them by
%     % 15% and drop labels on any axes not at the figure edge
%     figure;
%     ax=makesubplots(5,3,1:12);
%     ax=reshape(ax,3,4);
%     tax=ax';
%     axexpand(ax,15);
%     nolabels(tax(5:12),'y');
%     nolabels(ax(1:9),'x');
%     noticks(tax(5:12),'y');
%     noticks(ax(1:9),'x');
%     th=supertitle(ax,'This is a sooooooooooooooooooooooooper title!');
%     ax0=get(th,'parent'); % make title,colorbar,ylabel share same axis
%     superxlabel(ax0,'This is a sooooooooooooooooooooooooper xlabel!');
%     superylabel(ax0,'This is a sooooooooooooooooooooooooper ylabel!');
%     cb=supercolorbar(ax,'location','south');
%     cpos=get(cb,'position');
%     set(cb,'position',[cpos(1) cpos(2)-.15 cpos(3) cpos(4)/2]);
%     set(cb,'xaxislocation','bottom');
%
%    See also: NOTICKS, MAKESUBPLOTS, COMPACTAXES, AXEXPAND, AXMOVE,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR

%     Version History:
%        Aug.  4, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  4, 2010 at 12:25 GMT

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
    otherwise
        error('seizmo:nolabels:badInput',...
            'XY must be either ''x'' ''y'' or ''xy''!');
end

end
