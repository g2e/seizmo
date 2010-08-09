function notitles(ax)
%NOTITLES    Removes titles from specified axes
%
%    Usage:    notitles(ax)
%
%    Description:
%     NOTITLES(AX) removes titles from the specified axes.
%
%    Notes:
%
%    Examples:
%     % make and destroy a title
%     figure; axes;
%     title('hi...bye')
%     notitles
%
%    See also: NOLABELS, NOTICKS, NOCOLORBARS, SETFONTS, MAKESUBPLOTS,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              AXSTRETCH, AXEXPAND, AXMOVE

%     Version History:
%        Aug.  5, 2010 - initial version
%        Aug.  8, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% defaults
if(nargin<1 || isempty(ax)); ax=gca; end

% check axes
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:notitles:badInput',...
        'AX must be a valid axes handle!');
end

% remove titles
tmp=get(ax,'title');
if(iscell(tmp)); tmp=cat(1,tmp{:}); end
set(tmp,'string',[]);

end
