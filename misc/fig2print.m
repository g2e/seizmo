function fig2print(fh,orientation)
%FIG2PRINT    Adjusts figures to be as printed (aka print preview-ish)
%
%    Usage:    fig2print(fh,'orient')
%
%    Description:
%     FIG2PRINT(FH,'ORIENT') resizes the figure(s) given by the handle(s)
%     in FH so they match the printed size.  The 'ORIENT' input sets the
%     paper orientation.  Valid orientations are 'tall', 'portrait',
%     'landscape', 'ftall', 'fportrait', & 'flandscape'.  The last three
%     orientations do not have the quarter inch margins like the others.
%     Note that the figure itself is not necessarily the size of the print
%     (but the plots are the size they will be on the print).
%
%    Notes:
%
%    Examples:
%     % Setup a figure for landscape printing:
%     fig2print(gcf,'flandscape');
%
%    See also: PRINT, ORIENT, NOINVERT

%     Version History:
%        Jan.  5, 2002 - initial version
%        Aug.  4, 2010 - altered to match SEIZMO coding style
%        Apr.  4, 2012 - minor doc update
%
%     Written by Frederik J Simons (fjsimons-at-alum.mit.edu)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check figure handles exist
if(any(~ishandle(fh)) || any(~strcmpi('figure',get(fh,'type'))))
    error('seizmo:fig2print:badHandle',...
        'FH must be a valid figure handle!');
end

% check orientation
if(~ischar(orientation) || size(orientation,1)~=1 ...
        || ~isvector(orientation))
    error('seizmo:fig2print:badOrient',...
        'ORIENT must be a string!');
end

% loop over figures
for i=1:numel(fh)
    switch lower(orientation)
        case {'tall' 'portrait' 'landscape'}
            % just pass to Matlab's orient
            orient(fh(i),orientation);
        case {'ftall' 'fportrait'}
            %
            pu=get(fh(i),'PaperUnits');
            set(fh(i),'PaperUnits','normalized');
            set(fh(i),'Paperorientation','portrait',...
                'paperposition',[0 0 1 1]);
            set(fh(i),'PaperUnits',pu);
        case {'flandscape'}
            %
            pu=get(fh(i),'PaperUnits');
            set(fh(i),'PaperUnits','normalized');
            set(fh(i),'Paperorientation','landscape',...
                'paperposition',[0 0 1 1]);
            set(fh(i),'PaperUnits',pu);
        otherwise
            error('seizmo:fig2print:badOrient',...
                'Paper orientation unknown: %s',orientation);
    end
    
    % now expand figure to paper size
    % - get paper position
    % - get figure position in paper units
    % - expand figure
    ppos=get(fh(i),'PaperPosition');
    pu=get(fh(i),'PaperUnits');
    su=get(fh(i),'Units'); 
    set(fh(i),'Units',pu);
    spos=get(fh(i),'Position');
    set(fh(i),'Position',[spos(1)-(ppos(3)-spos(3)) ...
        spos(2)-(ppos(4)-spos(4)) ppos(3) ppos(4)]);
    set(fh(i),'Units',su);
end

end
