function []=showflags(ax,lgc)
%SHOWFLAGS    Toggle visibility of marker flags in SEIZMO plot
%
%    Usage:    showflags(ax)
%              showflags(ax,lgc)
%
%    Description:
%     SHOWFLAGS(AX) toggles the visibility of marker flags in the axes in
%     AX while ignoring the marker visibility.  The axes in AX must have
%     been drawn by either PLOT0, PLOT1, or RECORDSECTION.  If AX is not
%     given or is empty, SHOWFLAGS uses the current axes.
%
%     SHOWFLAGS(AX,LGC) forces the visibility to LGC.  This is useful for
%     synchronizing the visibility of all marker flags in AX.  LGC may be
%     TRUE/FALSE or 'ON'/'OFF'.
%
%    Notes:
%     - Calls DRAWMARKERS if it has not been called yet for the given axes.
%       This is done if there are no flags.
%
%    Examples:
%     % force all flags to be invisible (leaving the markers visible):
%     ax=plot0(data,'markers',true);
%     showflags(ax,false);
%
%    See also: DRAWMARKERS, PLOT0, PLOT1, RECORDSECTION, SHOWFLAGS,
%              STRETCHMARKERS, FLIPFLAGS

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
    error('seizmo:showflags:badInput',...
        'AXES must be an array of axes handles!');
end

% loop over axes
for i=1:numel(ax)
    % get flags
    hflag=findobj(ax(i),'-regexp','tag','flag$');
    
    % draw markers and flags if there are none
    if(isempty(hflag))
        [hmark,hflag]=drawmarkers(ax(i));
    end
    if(isempty(hflag)); continue; end
    
    % get current visibility if flag not given
    if(nargin<2 || isempty(lgc))
        % get visibility
        hfv=get(hflag,'visible');
        
        % get unique visibility
        [uhfv,ihfv,ihfv]=unique(hfv);
        
        % do common ones together
        for j=1:numel(uhfv)
            if(strcmpi(uhfv(j),'on'))
                set(hflag(ihfv==j),'visible','off');
            else
                set(hflag(ihfv==j),'visible','on');
            end
        end
    else
        if(islogical(lgc) && isscalar(lgc))
            % convert to proper string
            if(lgc)
                lgc='on';
            else
                lgc='off';
            end
            
            % set visibility
            set(hflag,'visible',lgc);
        elseif(isstring(lgc) && any(strcmpi(lgc,{'on' 'off'})))
            % set visibility
            set(hflag,'visible',lgc);
        else
            error('seizmo:showflags:badInput',...
                'LGC must be TRUE/FALSE or ''ON''/''OFF''!');
        end
    end
end

end
