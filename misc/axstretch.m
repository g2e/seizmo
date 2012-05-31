function axstretch(ax,varargin)
%AXSTRETCH    Stretches a set of axes as a group (resizing+moving)
%
%    Usage:    axstretch(ax,'direction',percent)
%              axstretch(ax,'*direction',percent)
%              axstretch(ax,'dir1',percent1,'dir2',percent2,...)
%
%    Description:
%     AXSTRETCH(AX,'DIRECTION',PERCENT) stretches a group of axes given
%     by the handles in AX towards the specified direction 'DIRECTION' by
%     PERCENT.  'DIRECTION' should be one of the following:
%      'left'   'l'  stretches to the left
%      'right'  'r'  stretches to the right
%      'up'     'u'  stretches upwards
%      'down'   'd'  stretches downwards
%      'horz'   'h'  stretches horizontally (preserving group center)
%      'vert'   'v'  stretches vertically (preserving group center)
%     The stretch is specified in percent of the total width/height of
%     the group of axes.  Positive PERCENT indicates expansion while
%     negative will contract the axes group.
%
%     AXSTRETCH(AX,'*DIRECTION',PERCENT) will move the axes in AX without
%     resizing them to satisfy the stretching requirements.
%
%     AXSTRETCH(AX,'DIR1',PERC1,'DIR2',PERC2,...) allows stretching in
%     multiple directions.  Note that the percentages are based on the
%     original axes.
%
%    Notes:
%
%    Examples:
%     % Make a group of subplots and stretch them downwards:
%     figure;
%     ax=makesubplots(5,3,1:12);
%     axstretch(ax,'down',20);
%
%    See also: AXMOVE, AXEXPAND, MAKESUBPLOTS, NOLABELS, NOTICKS, NOTITLES,
%              NOCOLORBARS, SUPERTITLE, SUPERXLABEL, SUPERYLABEL,
%              SUPERCOLORBAR

%     Version History:
%        Aug.  5, 2010 - initial version
%        Aug.  7, 2010 - added more options, flipped meaning of percent
%                        sign, changed name to axstretch
%        Aug.  8, 2010 - doc update
%        Dec. 17, 2010 - fixed buggy calls to min/max
%        May   4, 2012 - stretch called w/ multiple commands fixed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   4, 2012 at 12:25 GMT

% todo:
% - maybe add units triggers

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:axstretch:badPair',...
        'Incomplete Direction/Percentage pair!');
end

% check axes
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:axstretch:badInput',...
        'AX must be a valid axes handle!');
end

% check direction & percentage
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:axstretch:badInput',...
        'DIRECTION must be a string!');
elseif(any(~cellfun('isreal',varargin(2:2:end))) ...
        || any(cellfun('prodofsize',varargin(2:2:end))~=1) ...
        || any(isinf([varargin{2:2:end}])))
    error('seizmo:axstretch:badInput',...
        'PERCENT must be a finite valued scalar!');
end

% get position of "group axis"
lbwh=get(ax,'position');
if(iscell(lbwh)); lbwh=cat(1,lbwh{:}); end
rt=lbwh(:,1:2)+lbwh(:,3:4); % get rt for lbwh
gpos=[min(lbwh(:,1:2),[],1) max(rt,[],1)]; % LBRT
gpos(5:6)=gpos(3:4)-gpos(1:2); % LBRTWH

% loop over pairs
nax=numel(ax);
for i=1:2:nargin-1
    % percent to fraction
    frac=varargin{i+1}/100;
    
    switch lower(varargin{i})
        case {'left' 'l'}
            % compact/expand horizontally preserving rightmost edge
            newpos=[lbwh(:,1)+(lbwh(:,1)-gpos(3)).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1+frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'right' 'r'}
            % compact/expand horizontally preserving leftmost edge
            newpos=[lbwh(:,1)+(lbwh(:,1)-gpos(1)).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1+frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'up' 'u'}
            % compact/expand vectically preserving lowermost edge
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-gpos(2)).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1+frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'down' 'd'}
            % compact/expand vectically preserving uppermost edge
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-gpos(4)).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1+frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'horizontally' 'horz' 'h'}
            % compact/expand horizontally preserving center
            newpos=[lbwh(:,1)+(lbwh(:,1)-(gpos(1)+gpos(3))/2).*frac ...
                    lbwh(:,2) ...
                    lbwh(:,3).*(1+frac) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'vertically' 'vert' 'v'}
            % compact/expand vectically preserving center
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-(gpos(2)+gpos(4))/2).*frac ...
                    lbwh(:,3) ...
                    lbwh(:,4).*(1+frac)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*left' '*l'}
            % compact/expand horizontally preserving rightmost edge
            
            % weights
            lw=abs(1./((gpos(1)-lbwh(:,1))./gpos(5))); % 1 to inf
            rw=abs(1./((gpos(3)-rt(:,1))./gpos(5))); % 1 to inf
            
            % push factor
            push=pushfactor(lw,rw).*lbwh(:,3).*frac;
            
            % stretch left without axis expansion
            newpos=[lbwh(:,1)+(lbwh(:,1)-gpos(3)).*frac+push ...
                    lbwh(:,2) ...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*right' '*r'}
            % compact/expand horizontally preserving leftmost edge
            
            % weights
            lw=abs(1./((gpos(1)-lbwh(:,1))./gpos(5))); % 1 to inf
            rw=abs(1./((gpos(3)-rt(:,1))./gpos(5))); % 1 to inf
            
            % push factor
            push=pushfactor(lw,rw).*lbwh(:,3).*frac;
            
            % stretch right without axis expansion
            newpos=[lbwh(:,1)+(lbwh(:,1)-gpos(1)).*frac+push ...
                    lbwh(:,2) ...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*up' '*u'}
            % compact/expand vectically preserving lowermost edge
            
            % weights
            bw=abs(1./((gpos(2)-lbwh(:,2))./gpos(6))); % 1 to inf
            tw=abs(1./((gpos(4)-rt(:,2))./gpos(6))); % 1 to inf
            
            % push factor
            push=pushfactor(bw,tw).*lbwh(:,4).*frac;
            
            % stretch up without axis expansion
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-gpos(2)).*frac+push ...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*down' '*d'}
            % compact/expand vectically preserving uppermost edge
            
            % weights
            bw=abs(1./((gpos(2)-lbwh(:,2))./gpos(6))); % 1 to inf
            tw=abs(1./((gpos(4)-rt(:,2))./gpos(6))); % 1 to inf
            
            % push factor
            push=pushfactor(bw,tw).*lbwh(:,4).*frac;
            
            % stretch down without axis expansion
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-gpos(4)).*frac+push ...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*horizontally' '*horz' '*h'}
            % compact/expand horizontally preserving center
            
            % weights
            lw=abs(1./((gpos(1)-lbwh(:,1))./gpos(5))); % 1 to inf
            rw=abs(1./((gpos(3)-rt(:,1))./gpos(5))); % 1 to inf
            
            % push factor
            push=pushfactor(lw,rw).*lbwh(:,3).*frac;
            
            % stretch horizontally without axis expansion
            newpos=[lbwh(:,1)+(lbwh(:,1)-(gpos(1)+gpos(3))/2).*frac+push...
                    lbwh(:,2) ...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        case {'*vertically' '*vert' '*v'}
            % compact/expand vectically preserving center
            
            % weights
            bw=abs(1./((gpos(2)-lbwh(:,2))./gpos(6))); % 1 to inf
            tw=abs(1./((gpos(4)-rt(:,2))./gpos(6))); % 1 to inf
            
            % push factor
            push=pushfactor(bw,tw).*lbwh(:,4).*frac;
            
            % stretch down without axis expansion
            newpos=[lbwh(:,1) ...
                    lbwh(:,2)+(lbwh(:,2)-(gpos(2)+gpos(4))/2).*frac+push...
                    lbwh(:,3) ...
                    lbwh(:,4)];
            for j=1:nax; set(ax(j),'position',newpos(j,:)); end
        otherwise
            error('seizmo:axstretch:badInput',...
                'Unknown DIRECTION: %s',varargin{i});
    end
    
    % update lbwh, gpos, rt
    lbwh=newpos;
    rt=lbwh(:,1:2)+lbwh(:,3:4);
    gpos=[min(lbwh(:,1:2),[],1) max(rt,[],1)]; % LBRT
    gpos(5:6)=gpos(3:4)-gpos(1:2); % LBRTWH
end

end

function [push]=pushfactor(lw,rw)
% 0 means use none of the stretch factor to shift from l/b corner
% 1 means use all of the stretch factor to shift to r/t corner
% .5 will use half of the stretch factor to shift to center

% fix inf
linf=isinf(lw);
rinf=isinf(rw);
lw(linf & rinf)=1;
rw(linf & rinf)=1;
lw(linf & ~rinf)=1;
rw(linf & ~rinf)=0;
lw(~linf & rinf)=0;
rw(~linf & rinf)=1;

% norm
nw=lw+rw;
lw=lw./nw;
rw=rw./nw;

% get push factor
push=(1+rw-lw)./2; % -1 to 1

end
