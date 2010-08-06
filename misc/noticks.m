function noticks(ax,xy)
%NOTICKS    Removes ticks and tick labels from axes
%
%    Usage:    noticks()
%              noticks(ax)
%              noticks(ax,'x')
%              noticks(ax,'y')
%
%    Description:
%     NOTICKS() removes x/y ticks (and tick labels) from the current axis.
%
%     NOTICKS(AX) removes x/y ticks+labels from the specified axis/axes.
%
%     NOTICKS(AX,'X') only removes the x-axis ticks and their labels.
%
%     NOTICKS(AX,'Y') only removes the y-axis ticks and their labels.
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
%    See also: NOLABELS, MAKESUBPLOTS, COMPACTAXES, AXEXPAND, AXMOVE,
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
    error('seizmo:noticks:badInput',...
        'AX must be a valid axes handle!');
elseif(~ischar(xy))
    error('seizmo:noticks:badInput',...
        'XY must be either ''x'' ''y'' or ''xy''!');
end

% implement on axis of choice
switch lower(xy)
    case 'xy'
        set(ax,'xtick',[],'ytick',[]);
    case 'x'
        set(ax,'xtick',[]);
    case 'y'
        set(ax,'ytick',[]);
    otherwise
        error('seizmo:noticks:badInput',...
            'XY must be either ''x'' ''y'' or ''xy''!');
end

end
