function []=key2zoompan(key,ax)
%KEY2ZOOMPAN    Zooms or pans axis using keyboard
%
%    Usage:    key2zoompan(key)
%              key2zoompan(key,ax)
%
%    Description:
%     KEY2ZOOMPAN(KEY) will possibly zoom or pan the current axis depending
%     on the keycode given by KEY.  A keycode is an integer corresponding
%     to a specific key on the keyboard.  You can obtain a keycode via the
%     button variable returned by GINPUT.  Actions are currently:
%      KEY          KEYCODE    ACTION
%       up-arrow     30         Pan up 25%
%       down-arrow   31         Pan down 25%
%       left-arrow   28         Pan left 25%
%       right-arrow  29         Pan right 25%
%       i            105        Zoom in 2x
%       o            111        Zoom out 2x
%       x            120        Toggle x-axis only zoom on/off (def is off)
%       y            121        Toggle y-axis only zoom on/off (def is off)
%       r            114        Zoom out to default zoom level
%       s            115        Set default zoom level to current
%
%     KEY2ZOOMPAN(KEY,AX) uses the axis given by AX.  Note that if the axis
%     handle changes, the x/y-axis zoom states are lost for the previous
%     axis.
%
%    Notes:
%     - If the caller to KEY2ZOOMPAN changes then the x/y-axis zoom states
%       from the previous caller are lost.
%
%    Examples:
%     % Plot some data and zoom/pan the axis using the keyboard:
%     ax=scatter(rand(100,1),rand(100,1));
%     while(true)
%         [x,y,b]=ginput(1);
%         if(ishandle(ax) && ax==gca)
%              key2zoompan(b,ax);
%         end
%     end
%
%    See also: USERWINDOW, USERTAPER, USERSNR, USERWINNOW, SELECTRECORDS,
%              SELECTCLUSTERS, PICKSTACK, SELECTCLUSTERS_AND_GROUNDUNITS,
%              USERCLUSTER, ERRCUT, ARRCUT, AMPCUT, SNRCUT, POPCUT

%     Version History:
%        Jan.  6, 2011 - initial version
%        Jan.  7, 2011 - use axis handle for zoom calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  7, 2011 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default axis (current)
if(nargin<2 || isempty(ax)); ax=gca; end

% check inputs
if(~isscalar(key) || ~isreal(key) || key~=fix(key) || key<1)
    error('seizmo:key2zoompan:badInput',...
        'KEY must be a integer-valued scalar!');
elseif(~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    error('seizmo:key2zoompan:badInput',...
        'AX must be a valid axes handle!');
end

% whos calling
caller=star69;

% save & default x/y axis zoom states
persistent xzoom yzoom oldax oldcaller
if(isempty(xzoom) || ax~=oldax || ~strcmp(caller,oldcaller))
    xzoom=false;
end
if(isempty(yzoom) || ax~=oldax || ~strcmp(caller,oldcaller));
    yzoom=false;
end

% action based on key
switch key
    case 105 % i
        if(xzoom)
            zoom(ax,'xon');
            zoom(ax,2);
            zoom(ax,'on'); zoom(ax,'off');
        elseif(yzoom)
            zoom(ax,'yon');
            zoom(ax,2);
            zoom(ax,'on'); zoom(ax,'off');
        else
            zoom(ax,2);
        end
    case 111 % o
        if(xzoom)
            zoom(ax,'xon');
            zoom(ax,.5);
            zoom(ax,'on'); zoom(ax,'off');
        elseif(yzoom)
            zoom(ax,'yon');
            zoom(ax,.5);
            zoom(ax,'on'); zoom(ax,'off');
        else
            zoom(ax,.5);
        end
    case 120 % x
        if(xzoom)
            xzoom=false;
        else
            xzoom=true;
            yzoom=false;
        end
    case 121 % y
        if(yzoom)
            yzoom=false;
        else
            yzoom=true;
            xzoom=false;
        end
    case 115 % s
        zoom(ax,'reset');
    case 114 % r
        zoom(ax,'out');
    case 30 % up arrow
        ylims=ylim(ax);
        switch get(ax,'ydir')
            case 'reverse'
                ylim(ax,ylims-(ylims(2)-ylims(1))/4);
            otherwise
                ylim(ax,ylims+(ylims(2)-ylims(1))/4);
        end
    case 31 % down arrow
        ylims=ylim(ax);
        switch get(ax,'ydir')
            case 'reverse'
                ylim(ax,ylims+(ylims(2)-ylims(1))/4);
            otherwise
                ylim(ax,ylims-(ylims(2)-ylims(1))/4);
        end
    case 29 % right arrow
        xlims=xlim(ax);
        xlim(ax,xlims+(xlims(2)-xlims(1))/4);
    case 28 % left arrow
        xlims=xlim(ax);
        xlim(ax,xlims-(xlims(2)-xlims(1))/4);
    otherwise
        % nothing
end

% set oldax & oldcaller
oldax=ax;
oldcaller=caller;

end
