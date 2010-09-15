function []=stretchmarkers(ax,perc)
%STRETCHMARKERS    Stretches markers vertically by specified percent
%
%    Usage:    stretchmarkers(ax,percent)
%
%    Description:
%     STRETCHMARKERS(AX,PERCENT) adjusts the vertical extent of markers in
%     axes AX by the specified percent PERCENT.  The axes in AX must have
%     been drawn by either PLOT0, PLOT1, or RECORDSECTION.  PERCENT must be
%     a real-valued scalar.  A PERCENT of 50 will increase the marker's
%     height by 50%.
%
%    Notes:
%
%    Examples:
%     % reduce the extent of markers by 50%
%     plot0(data,'markers',true);
%     stretchmarkers([],-50);
%
%    See also: DRAWMARKERS, SHOWMARKERS, SHOWFLAGS, FLIPFLAGS

%     Version History:
%        Sep. 14, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2010 at 23:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% use current axes if none given
if(isempty(ax)); ax=gca; end

% check all inputs are axes
if(isempty(ax) || ~isreal(ax) || any(~ishandle(ax)) ...
        || any(~strcmp('axes',get(ax,'type'))))
    error('seizmo:stretchmarkers:badInput',...
        'AXES must be an array of axes handles!');
end

% check percent is valid
if(~isreal(perc) || ~isscalar(perc))
    error('seizmo:stretchmarkers:badInput',...
        'PERCENT must be a real-valued scalar!');
end

% loop over axes
for i=1:numel(ax)
    % get markers and flags
    hmark=findobj(ax(i),'-regexp','tag','marker$');
    hflag=findobj(ax(i),'-regexp','tag','flag$');
    
    % skip if no markers left
    if(isempty(hmark)); continue; end
    
    % get top/bottom of marker
    ydata=get(hmark,'ydata');
    if(~iscell(ydata)); ydata={ydata}; end
    
    % loop over each marker
    for j=1:numel(ydata)
        % adjust ydata
        ydata{j}=sum(ydata{j})/2+[-1 1].*diff(ydata{j}).*(100+perc)./200;
        set(hmark(j),'ydata',ydata{j});
    end
    
    % loop over each flag
    for j=1:numel(hflag)
        tmp=get(hflag(j),'userdata');
        if(~ishandle(tmp.marker)); continue; end
        xdata=get(tmp.marker,'xdata');
        ydata=get(tmp.marker,'ydata');
        set(hflag(j),...
            'position',[xdata(1) ydata(1)+diff(ydata)*tmp.mast/100]);
    end
end

end
