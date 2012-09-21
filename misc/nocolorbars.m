function nocolorbars(ax)
%NOCOLORBARS    Removes colorbars associated with the specified axes
%
%    Usage:    nocolorbars(ax)
%
%    Description:
%     NOCOLORBARS(AX) removes colorbars from the specified axes.
%
%    Notes:
%     - pass 'align' when making a subplot to avoid poor colorbar behavior
%
%    Examples:
%     % make and destroy a colorbar
%     figure; axes; colorbar;
%     nocolorbars;
%
%    See also: NOLABELS, NOTICKS, NOTITLES, SETFONTS, MAKESUBPLOTS,
%              SUPERTITLE, SUPERXLABEL, SUPERYLABEL, SUPERCOLORBAR,
%              AXSTRETCH, AXEXPAND, AXMOVE

%     Version History:
%        Aug.  5, 2010 - initial version
%        Aug.  8, 2010 - doc update
%        Oct. 11, 2010 - delete colorbar handles
%        Sep.  5, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2012 at 12:25 GMT

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

% remove colorbars
for i=1:numel(ax)
    cbh=findobj(get(ax(i),'parent'),'tag','Colorbar','axes',ax(i));
    if(~isempty(cbh)); delete(cbh); end
end

end
