function axexpand(ax,wf,hf)
%AXEXPAND    Expands axes by factor
%
%    Usage:    axexpand()
%              axexpand(ax)
%              axexpand(ax,factor)
%              axexpand(ax,widthfac,heightfac)
%
%    Description:
%     AXEXPAND() expands the current axis by 10%. This just alters the axis
%     'position' property.  Both the vertical and horizontal sizes are
%     altered.
%
%     AXEXPAND(AX) expands the specified axes in AX by 10%.
%
%     AXEXPAND(AX,FACTOR) specifies the expansion factor in percent.  This
%     is applied to both the width and height.  0% is no change.  Negative
%     values will cause shrinkage.  May be a vector if AX is a vector.
%
%     AXEXPAND(AX,WIDTHFAC,HEIGHTFAC) sets separate factors for the width
%     and height of an axes.  Values are in percent where 0% is no change.
%
%    Notes:
%     - To undo an expansion: NEWFAC=(1/(1+OLDFAC/100)-1)*100
%
%    Examples:
%     % Make a figure with 4x3 arrangement of subplots, then expand them by
%     % 15% and drop labels on any axes not at the figure edge:
%     figure;
%     ax=makesubplots(5,3,1:12);
%     ax=reshape(ax,3,4);
%     tax=ax';
%     axexpand(ax,15);
%     nolabels(tax(5:12),'y');
%     nolabels(ax(1:9),'x');
%
%    See also: AXSTRETCH, AXMOVE, MAKESUBPLOTS, NOLABELS, NOTICKS,
%              NOTITLES, NOCOLORBARS, SUPERTITLE, SUPERXLABEL, SUPERYLABEL,
%              SUPERCOLORBAR

%     Version History:
%        Aug.  4, 2010 - initial version
%        Aug.  8, 2010 - doc update
%        May   4, 2012 - minor doc fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   4, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% defaults
if(nargin<1 || isempty(ax)); ax=gca; end
if(nargin<2 || isempty(wf)); wf=10; end
if(nargin<3 || isempty(hf)); hf=wf; end

% check inputs
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:axexpand:badInput',...
        'AX must be a valid axes handle!');
elseif(nargin<3 && ~isreal(wf))
    error('seizmo:axexpand:badInput',...
        'FACTOR must be a real value in percent!');
elseif(nargin>2 && ~isreal(wf))
    error('seizmo:axexpand:badInput',...
        'WIDTHFAC must be a real value in percent!');
elseif(~isreal(hf))
    error('seizmo:axexpand:badInput',...
        'HEIGHTFAC must be a real value in percent!');
elseif(~isequalsizeorscalar(ax,wf,hf))
    error('seizmo:axexpand:badInput',...
        'All inputs must be equal sized or scalar!');
end

% expand inputs
[ax,wf,hf]=expandscalars(ax,wf,hf);

% expand axes
for i=1:numel(ax)
    % current position
    pos=get(ax(i),'position');
    
    % midpoint
    midx=pos(1)+pos(3)/2;
    midy=pos(2)+pos(4)/2;
    
    % new position
    neww=pos(3)*(1+wf(i)/100);
    newh=pos(4)*(1+hf(i)/100);
    set(ax(i),'position',[midx-neww/2 midy-newh/2 neww newh]);
end

end
