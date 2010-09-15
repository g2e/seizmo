function []=flipflags(ax,dir)
%FLIPFLAGS    Flips the flag direction in SEIZMO plots
%
%    Usage:    flipflags(ax)
%              flipflags(ax,dir)
%
%    Description:
%     FLIPFLAGS(AX) flips the direction of marker flags in the axes given
%     by AX.  The axes in AX must have been drawn by PLOT0, PLOT1 or
%     RECORDSECTION.  If AX is not given or is empty, FLIPFLAGS uses the
%     current axes.
%
%     FLIPFLAGS(AX,DIR) forces the flag direction to DIR.  This is useful
%     for synchronizing flag directions.  DIR may be 'LEFT', 'CENTER', or
%     'RIGHT' and indicates the flag position relative to the marker.
%
%    Notes:
%
%    Examples:
%     % Flip the flags to the left of the markers
%     ax=plot0(data,'markers',true);
%     flipflags(ax,'left');
%
%    See also: SHOWFLAGS, MASTFLAGS, SHOWMARKERS, DRAWMARKERS,
%              STRETCHMARKERS

%     Version History:
%        Sep. 14, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2010 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% use current axes if none given
if(nargin<1 || isempty(ax)); ax=gca; end

% check all inputs are axes
if(isempty(ax) || ~isreal(ax) || any(~ishandle(ax)) ...
        || any(~strcmp('axes',get(ax,'type'))))
    error('seizmo:showmarkers:badInput',...
        'AXES must be an array of axes handles!');
end

% loop over axes
for i=1:numel(ax)
    % get flags
    hflag=findobj(ax(i),'-regexp','tag','flag$');
    
    % skip if none
    if(isempty(hflag)); continue; end
    
    % get current direction if not given
    if(nargin<2 || isempty(dir))
        % get direction
        dir=get(hflag,'horizontalalignment');
        
        % get unique directions
        [udir,idir,idir]=unique(dir);
        
        % do common ones together
        for j=1:numel(udir)
            if(strcmpi(udir(j),'left'))
                set(hflag(idir==j),'horizontalalignment','right');
            elseif(strcmpi(udir(j),'right'))
                set(hflag(idir==j),'horizontalalignment','left');
            end
        end
    else
        if(isstring(dir) && any(strcmpi(dir,{'left' 'right' 'center'})))
            % switch left and right b/c we are indicating the which side
            % of the marker the flags are on, not which side of the flags
            % are touching the marker
            if(strcmpi(dir,'left'))
                dir='right';
            elseif(strcmpi(dir,'right'))
                dir='left';
            end
            set(hflag,'horizontalalignment',dir);
        else
            error('seizmo:flipflags:badInput',...
                'DIR must be ''left'',''center'', or ''right''!');
        end
    end
end

end
