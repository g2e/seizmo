function compactaxes(ax,varargin)
%COMPACTAXES    Compacts (resizes) a set of axes in a specific direction
%
%    Usage:    compactaxes(ax,'direction',percent)
%              compactaxes(ax,'dir1',percent1,'dir2',percent2,...)
%
%    Description:
%     COMPACTAXES(AX,'DIRECTION',PERCENT) compacts a group of axes given by
%     the handles in AX towards the specified direction 'DIRECTION' by
%     PERCENT.  'DIRECTION' should be one of the following:
%      'left'   'l'  compacts to the left
%      'right'  'r'  compacts to the right
%      'up'     'u'  compacts upwards
%      'down'   'd'  compacts downwards
%      'horz'   'h'  compacts horizontally (preserving group center)
%      'vert'   'v'  compacts vertically (preserving group center)
%     The compaction is specified in percent of the total width/height of
%     the group of axes.  Positive PERCENT indicates compaction while
%     negative will expand the axes group.
%
%     COMPACTAXES(AX,'DIR1',PERC1,'DIR2',PERC2,...) allows compaction in
%     multiple directions.  Note that the percentages are based on the
%     original axes.
%
%    Notes:
%
%    Examples:
%     % make a group of subplots and expand them downwards:
%      figure;
%      ax=makesubplots(5,3,1:12);
%      compactaxes(ax,'up',-20);
%
%    See also: AXMOVE, AXEXPAND

%     Version History:
%        Aug.  5, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  5, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:compactaxes:badPair',...
        'Incomplete Direction/Percentage pair!');
end

% check axes
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:compactaxes:badInput',...
        'AX must be a valid axes handle!');
end

% check direction & percentage
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:compactaxes:badInput',...
        'DIRECTION must be a string!');
elseif(any(~cellfun('isreal',varargin(2:2:end))) ...
        || any(cellfun('prodofsize',varargin(2:2:end))~=1))
    error('seizmo:compactaxes:badInput',...
        'PERCENT must be a scalar!');
end

% get position of "group axis"
lbwh=get(ax,'position');
if(iscell(lbwh)); lbwh=cat(1,lbwh{:}); end
gpos=[min(lbwh(:,1:2)) max(lbwh(:,1:2)+lbwh(:,3:4))]; % LBRT
gpos(5:6)=gpos(3:4)-gpos(1:2); % LBRTWH

% loop over pairs
nax=numel(ax);
for i=1:2:nargin-1
    % percent to fraction
    frac=varargin{i+1}/100;
    
    switch lower(varargin{i})
        case {'left' 'l'}
            % compact/expand horizontally preserving leftmost edge
            newpos=[lbwh(:,1)-(lbwh(:,1)-gpos(1)).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1-frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'right' 'r'}
            % compact/expand horizontally preserving rightmost edge
            newpos=[lbwh(:,1)-(lbwh(:,1)-gpos(3)).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1-frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'up' 'u'}
            % compact/expand vectically preserving uppermost edge
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)-(lbwh(:,2)-gpos(4)).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1-frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'down' 'd'}
            % compact/expand vectically preserving lowermost edge
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)-(lbwh(:,2)-gpos(2)).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1-frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'horizontally' 'horz' 'h'}
            % compact/expand horizontally preserving center
            newpos=[lbwh(:,1)-(lbwh(:,1)-(gpos(1)+gpos(3))/2).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1-frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'vertically' 'vert' 'v'}
            % compact/expand vectically preserving center
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)-(lbwh(:,2)-(gpos(2)+gpos(4))/2).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1-frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        otherwise
            error('seizmo:compactaxes:badInput',...
                'Unknown DIRECTION: %s',varargin{i});
    end
end

end
