function []=showmarkers(ax,lgc)
%SHOWMARKERS    Toggle visibility of markers & flags in SEIZMO plot
%
%    Usage:    showmarkers(ax)
%              showmarkers(ax,lgc)
%
%    Description:
%     SHOWMARKERS(AX) toggles the visibility of markers & marker flags in
%     the axes in AX.  The axes in AX must have been drawn by either PLOT0,
%     PLOT1, or RECORDSECTION.  If AX is not given or is empty, SHOWMARKERS
%     uses the current axes.
%
%     SHOWMARKERS(AX,LGC) forces the visibility to LGC.  This is useful for
%     synchronizing the visibility of all flags and markers in AX.  LGC may
%     be TRUE/FALSE or 'ON'/'OFF'.
%
%    Notes:
%     - Calls DRAWMARKERS if it has not been called yet for the given axes.
%       This is done if there are no markers or no flags.
%
%    Examples:
%     % force all flags and markers to be invisible
%     ax=plot0(data,'markers',true);
%     showmarkers(ax,false);
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
    error('seizmo:showmarkers:badInput',...
        'AXES must be an array of axes handles!');
end

% loop over axes
for i=1:numel(ax)
    % get markers and flags
    hmark=findobj(ax(i),'-regexp','tag','marker$');
    hflag=findobj(ax(i),'-regexp','tag','flag$');
    
    % draw markers and flags if not there
    if(isempty(hmark) || isempty(hflag))
        [hmark,hflag]=drawmarkers(ax(i));
    end
    if(isempty(hmark) || isempty(hflag)); continue; end
    
    % get current visibility if flag not given
    if(nargin<2 || isempty(lgc))
        % get visibility
        hmv=get(hmark,'visible');
        hfv=get(hflag,'visible');
        
        % get unique visibility
        [uhmv,ihmv,ihmv]=unique(hmv);
        [uhfv,ihfv,ihfv]=unique(hfv);
        
        % do common ones together
        for j=1:numel(uhmv)
            if(strcmpi(uhmv(j),'on'))
                set(hmark(ihmv==j),'visible','off');
            else
                set(hmark(ihmv==j),'visible','on');
            end
        end
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
            set(hmark,'visible',lgc);
            set(hflag,'visible',lgc);
        elseif(isstring(lgc) && any(strcmpi(lgc,{'on' 'off'})))
            % set visibility
            set(hmark,'visible',lgc);
            set(hflag,'visible',lgc);
        else
            error('seizmo:showmarkers:badInput',...
                'LGC must be TRUE/FALSE or ''ON''/''OFF''!');
        end
    end
end

end
