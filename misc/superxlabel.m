function [h]=superxlabel(ax,varargin)
%SUPERXLABEL    Makes an x-axis label spanning multiple axes
%
%    Usage:    superxlabel(ax,'text')
%              superxlabel(ax,'text','prop','val',...)
%              h=superxlabel(...)
%
%    Description:
%     SUPERXLABEL(AX,'TEXT') adds a label to the bottom of a group of axes
%     given in AX.  All axes in AX must share the same figure!
%
%     SUPERXLABEL(AX,'TEXT','PROP','VAL',...) sets the values for the
%     specified label properties.
%
%     H=SUPERXLABEL(...) returns the handle to the text object used as the
%     label.  To get handle to the invisible axis for the label use
%     GET(H,'PARENT').
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
%    See also: SUPERTITLE, SUPERYLABEL, SUPERCOLORBAR

%     Version History:
%        Aug.  5, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  5, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check inputs
if(isempty(ax) || (isreal(ax) && isscalar(ax)) ...
        || ischar(ax) || iscellstr(ax))
    % single or no axis - use xlabel
    h=xlabel(ax,varargin{:});
    set(h,'visible','on');
    return;
elseif(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))) ...
        || ~isscalar(unique(cell2mat(get(ax(:),'parent')))))
    error('seizmo:superxlabel:badInput',...
        'AX must be valid axes handles all in the same figure!');
end

% get position of new axis
lbwh=get(ax,'position');
if(iscell(lbwh)); lbwh=cat(1,lbwh{:}); end
newpos=[min(lbwh(:,1:2)) max(lbwh(:,1:2)+lbwh(:,3:4))]; % LBRT
newpos(3:4)=newpos(3:4)-newpos(1:2); % LBWH

% create axis, set xlabel, make invisible
ax=axes('position',newpos);
h=xlabel(ax,varargin{:});
set(ax,'visible','off');
set(h,'visible','on');

end
