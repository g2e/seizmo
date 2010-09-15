function []=mastflags(ax,mast)
%MASTFLAGS    Masts marker flags in SEIZMO plots
%
%    Usage:    mastflags(ax)
%              mastflags(ax,mast)
%
%    Description:
%     MASTFLAGS(AX) sets the marker flags in axes given by AX to half-mast.
%     The axes in AX must have been drawn by PLOT0, PLOT1 or
%     RECORDSECTION.  If AX is not given or is empty, MASTFLAGS uses the
%     current axes.
%
%     MASTFLAGS(AX,MAST) sets the flag mast to a custom position given by
%     MAST.  MAST must be a scalar from 0 to 100.  The default is 50.
%
%    Notes:
%
%    Examples:
%     % Set flags to the bottom of the marker
%     ax=plot0(data,'markervertalign','bottom');
%     mastflags(ax,0);
%
%    See also: SHOWFLAGS, FLIPFLAGS, SHOWMARKERS, DRAWMARKERS,
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
if(nargin<2 || isempty(mast)); mast=50; end

% check all inputs are axes
if(isempty(ax) || ~isreal(ax) || any(~ishandle(ax)) ...
        || any(~strcmp('axes',get(ax,'type'))))
    error('seizmo:showmarkers:badInput',...
        'AXES must be an array of axes handles!');
end

% check mast
if(~isreal(mast) || ~isscalar(mast))
    error('seizmo:showmarkers:badInput',...
        'MAST must be a real-valued scalar from 0 to 100!');
end

% loop over axes
for i=1:numel(ax)
    % get flags
    hmark=findobj(ax(i),'-regexp','tag','marker$');
    hflag=findobj(ax(i),'-regexp','tag','flag$');
    
    % skip if no markers/flags left
    if(isempty(hflag) || isempty(hmark)); continue; end
    
    % grab function type
    userdata=get(ax(i),'userdata');
    if(isstruct(userdata) ...
            && all(isfield(userdata,{'markers' 'function'})))
        % handle plot0 ydir being opposite
        if(strcmpi(userdata.function,'plot0'))
            mast=100-mast;
        end
    end
    
    % loop over each flag
    for j=1:numel(hflag)
        tmp=get(hflag(j),'userdata');
        tmp.mast=mast;
        if(~ishandle(tmp.marker))
            set(hflag(j),'userdata',tmp);
            continue;
        end
        xdata=get(tmp.marker,'xdata');
        ydata=get(tmp.marker,'ydata');
        set(hflag(j),...
            'position',[xdata(1) ydata(1)+diff(ydata)*mast/100],...
            'userdata',tmp);
    end
end

end
