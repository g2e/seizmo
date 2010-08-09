function setfonts(ax,varargin)
%SETFONTS    Sets font props for all text objects in the specified axes
%
%    Usage:    setfonts(ax,'prop1',val1,'prop2',val2,...)
%
%    Description:
%     SETFONTS(AX,'PROP1',VAL1,'PROP2',VAL2,...) passes property/value
%     pairs to all the text objects in the axes handles specified by AX.
%     This is useful for mass font changes & synchronization.
%
%    Notes:
%
%    Examples:
%     % change font size and weight for all axes in a figure
%      setfonts(get(gcf,'children'),'fontsize',16,'fontweight','bold')
%
%    See also: SET, FINDALL, NOTICKS, NOCOLORBARS, NOTITLES, NOLABELS,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              AXSTRETCH, AXEXPAND, AXMOVE, MAKESUBPLOTS

%     Version History:
%        Aug.  5, 2010 - initial version
%        Aug.  8, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:setfonts:badPair',...
        'Incomplete property/value pair!');
end

% check axes
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:setfonts:badInput',...
        'AX must be a valid axes handle!');
end

% first adjust axes themselves
set(ax,varargin{:});

% find all text objects in these axes
% and adjust text properties
set(findall(ax,'type','text'),varargin{:});

end
