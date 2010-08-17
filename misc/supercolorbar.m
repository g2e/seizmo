function [h]=supercolorbar(ax,varargin)
%SUPERCOLORBAR    Makes a colorbar spanning multiple axes
%
%    Usage:    supercolorbar(ax)
%              supercolorbar(ax,'prop','val',...)
%              h=supercolorbar(...)
%
%    Description:
%     SUPERCOLORBAR(AX) adds a colorbar to the right of a group of axes
%     given in AX.  All axes in AX must share the same figure and have the
%     same clim range!
%
%     SUPERCOLORBAR(AX,'PROP','VAL',...) sets the values for the colorbar
%     properties.  Using the 'PEER' property will override the peer axis
%     set here.
%
%     H=SUPERCOLORBAR(...) returns the handle to the colorbar axis.  To get
%     handle to the invisible axis for the colorbar use GET(H,'PARENT').
%
%    Notes:
%     - Note that this does not move the given axes to make space for the
%       colorbar.  In fact, we also make sure it doesn't move the invisible
%       axis it is peered with so things stay aligned.
%
%    Examples:
%     % make a figure with 4 2x2 groups of subplots and add super
%     % labeling and super colorbars to each group
%     fh=figure;
%     set(fh,'position',get(fh,'position').*[1 1 1.5 1.5]);
%     ax=makesubplots(5,5,submat(lind(5),1:2,[1 2 4 5]),'parent',fh);
%     ax=mat2cell(reshape(ax,4,4),[2 2],[2 2]);
%     for i=1:4
%         supertitle(ax{i},['super title ' num2str(i)]);
%         superxlabel(ax{i},['super xlabel ' num2str(i)]);
%         superylabel(ax{i},['super ylabel ' num2str(i)]);
%         supercolorbar(ax{i},'location','eastoutside');
%     end
%
%    See also: SUPERTITLE, SUPERXLABEL, SUPERYLABEL, MAKESUBPLOTS,
%              NOLABELS, NOTICKS, NOTITLES, NOCOLORBARS, AXMOVE, AXEXPAND,
%              AXSTRETCH

%     Version History:
%        Aug.  5, 2010 - initial version
%        Aug.  8, 2010 - move super axis below, tag & userdata used to
%                        replace on subsequent calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check inputs
if(isempty(ax) || ischar(ax))
    % no axis - use colorbar
    h=colorbar(varargin{:});
    set(h,'visible','on');
    return;
elseif(isreal(ax) && isscalar(ax))
    % single axis - use colorbar
    h=colorbar('peer',ax,varargin{:});
    set(h,'visible','on');
    return;
elseif(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))) ...
        || ~isscalar(unique(cell2mat(get(ax(:),'parent')))) ...
        || ~isvector(unique(cell2mat(get(ax(:),'clim')),'rows')))
    error('seizmo:supercolorbar:badInput',...
        ['AX must be valid axes handles all in the same figure ' ...
        '& sharing the same clim!']);
end
p=get(ax(1),'parent');

% get position of new axis
lbwh=get(ax,'position');
if(iscell(lbwh)); lbwh=cat(1,lbwh{:}); end
newpos=[min(lbwh(:,1:2)) max(lbwh(:,1:2)+lbwh(:,3:4))]; % LBRT
newpos(3:4)=newpos(3:4)-newpos(1:2); % LBWH

% create axis, set colorbar, move below & make invisible
sax=findobj(p,'type','axes','tag','super','userdata',ax);
if(isempty(sax))
    sax=axes('position',newpos);
    set(sax,'parent',p);
    pos=get(sax,'position'); % to beat colorbar
    set(sax,'clim',get(ax(1),'clim'));
    h=colorbar('peer',sax,varargin{:});
    kids=get(p,'children');
    set(p,'children',[kids(2:end); kids(1)]);
    set(sax,'visible','off','tag','super','userdata',ax);
    set(sax,'position',pos);
    set(h,'visible','on');
else
    pos=get(sax,'position'); % to beat colorbar
    set(sax,'clim',get(ax(1),'clim'));
    h=colorbar('peer',sax,varargin{:});
    set(sax,'position',pos);
    set(h,'visible','on');
end

end
