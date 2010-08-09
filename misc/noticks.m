function noticks(ax,xy)
%NOTICKS    Removes ticks and tick labels from axes
%
%    Usage:    noticks()
%              noticks(ax)
%              noticks(ax,'x')
%              noticks(ax,'y')
%
%    Description:
%     NOTICKS() removes x/y ticks and tick labels from the current axis.
%
%     NOTICKS(AX) removes x/y ticks and tick labels from the specified
%     axis or axes.
%
%     NOTICKS(AX,'X') only removes the x-axis ticks and their labels.
%
%     NOTICKS(AX,'Y') only removes the y-axis ticks and their labels.
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
%     nolabels(tax(5:12),'y');
%     nolabels(ax(1:9),'x');
%     noticks(tax(5:12),'y');
%     noticks(ax(1:9),'x');
%
%    See also: NOLABELS, NOCOLORBARS, NOTITLES, SETFONTS, MAKESUBPLOTS
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              AXSTRETCH, AXEXPAND, AXMOVE

%     Version History:
%        Aug.  4, 2010 - initial version
%        Aug.  8, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 12:25 GMT

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
