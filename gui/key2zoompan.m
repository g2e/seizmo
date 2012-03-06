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
%      KEY           KEYCODE    ACTION
%       up-arrow      30         Pan up 25%
%       down-arrow    31         Pan down 25%
%       left-arrow    28         Pan left 25%
%       right-arrow   29         Pan right 25%
%       i             105        Zoom in 2x
%       o             111        Zoom out 2x
%       x             120        Toggle x-axis only zoom on/off (def=off)
%       y             121        Toggle y-axis only zoom on/off (def=off)
%       r             114        Zoom out to default zoom level
%       s             115        Set default zoom level to current
%       h             104        Help window showing allowed keys
%
%     KEY2ZOOMPAN(KEY,AX) uses the axis given by AX.
%
%    Notes:
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
%        Mar.  6, 2012 - added h key for help window, now works while
%                        switching back and forth between axes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2012 at 11:00 GMT

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

% get states
persistent oldax oldzoom h
if(any(ax==oldax))
    xzoom=oldzoom(ax==oldax,1);
    yzoom=oldzoom(ax==oldax,2);
else
    xzoom=false;
    yzoom=false;
end

% action based on key
switch key
    case 104 % h
        k2zphelp={
            'KEY                ACTION'
            'up-arrow           Pan up 25%'
            'down-arrow         Pan down 25%'
            'left-arrow         Pan left 25%'
            'right-arrow        Pan right 25%'
            'i                  Zoom in 2x'
            'o                  Zoom out 2x'
            'x                  Toggle x-axis only zoom on/off (def=off)'
            'y                  Toggle y-axis only zoom on/off (def=off)'
            'r                  Zoom out to default zoom level'
            's                  Set default zoom level to current'
            'h                  This help window'};
        
        % does a help dialog window already exist?
        if(ishandle(h))
             % bring to front
            figure(h);
        else % new help dialog window
            h=helpdlg(k2zphelp,'KEY2ZOOMPAN CONTROL');
        end
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

% set oldax & oldzoom
if(any(ax==oldax))
    oldzoom(ax==oldax,:)=[xzoom yzoom];
else
    oldax=[oldax; ax];
    oldzoom=[oldzoom; xzoom yzoom];
end

% clear out dead handles
good=ishandle(oldax);
oldax=oldax(good);
oldzoom=oldzoom(good,:);

end
