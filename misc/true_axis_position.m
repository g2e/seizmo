function [pos]=true_axis_position(ax)
%TRUE_AXIS_POSITION    Returns plotted axis position (eg for "axis tight")
%
%    Usage:    pos=true_axis_position
%              pos=true_axis_position(ax)
%
%    Description:
%     POS=TRUE_AXIS_POSITION returns the true positioning of the current
%     axis which may differ from the axis positioning from GET, depending
%     on axis limits, data aspect ratio, and plot box aspect ratio.  The
%     position is returned in the same units as the those used to define
%     the axis.  This function can only be used for a 2D plot.  POS is a
%     4-element position vector like returned by: get(gca,'position').
%
%     POS=TRUE_AXIS_POSITION(AX) returns the position of the axis specified
%     by AX.
%
%    Notes:
%
%    Examples:
%     % Plot a line, adjust and show the difference between
%     % what this function returns & the position property:
%     plot([1 2],[1 10])
%     axis equal
%     get(gca,'position')
%     true_axis_position
%
%    See also: GET, SET, AXIS, AXES

%     Version History:
%        Jan. 13, 2006 - initial version
%        Feb. 14, 2006 - Clarified use of another File Exchange function
%                        (getInUnits.m).
%        May  24, 2007 - Fixed incorrect behavior when aspect ratios set to
%                        'auto', and some bugs regarding plot box aspect
%                        ratio.
%        May  26, 2010 - Small rewrite to remove dependency on external
%                        functions
%        Oct. 19, 2011 - Now defaults to current axis if no input provided.
%                        Also updated example image to a simpler one.
%        May.  3, 2012 - documentation changes, added to seizmo/misc
%        Oct. 15, 2012 - minor doc update
%
%     Written by Kelly Kearney (FEX ID: 9615)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2012 at 18:35 GMT

% todo:

% check input
error(nargchk(0,1,nargin));
if(nargin<1); ax=gca; end

% check handle
if(~isscalar(ax))
    error('seizmo:true_axis_position:badInput',...
        'AX must be a scalar axis handle!');
elseif(~ishghandle(ax,'axes'))
    error('seizmo:true_axis_position:badInput',...
        'AX must be an axis handle!');
end

% get position of axis in pixels
currunit=get(ax,'units');
set(ax,'units','pixels');
axispos=get(ax,'position');
set(ax,'units',currunit);

% calculate box position based axis limits and aspect ratios
darismanual=strcmpi(get(ax,'dataaspectratiomode'),'manual');
pbarismanual=strcmpi(get(ax,'plotboxaspectratiomode'),'manual');
if(~darismanual && ~pbarismanual)
    % everything is automatic so position property is right
    pos=axispos;
else % need to account for aspect ratios
    dx=diff(get(ax,'xlim'));
    dy=diff(get(ax,'ylim'));
    dar=get(ax,'dataaspectratio');
    pbar=get(ax,'plotboxaspectratio');

    limdarratio=(dx/dar(1))/(dy/dar(2));
    pbarratio=pbar(1)/pbar(2);
    axisratio=axispos(3)/axispos(4);

    if(darismanual)
        if(limdarratio>axisratio)
            pos(1)=axispos(1);
            pos(3)=axispos(3);
            pos(4)=axispos(3)/limdarratio;
            pos(2)=(axispos(4)-pos(4))/2+axispos(2);
        else
            pos(2)=axispos(2);
            pos(4)=axispos(4);
            pos(3)=axispos(4)*limdarratio;
            pos(1)=(axispos(3)-pos(3))/2+axispos(1);
        end
    elseif(pbarismanual)
        if(pbarratio>axisratio)
            pos(1)=axispos(1);
            pos(3)=axispos(3);
            pos(4)=axispos(3)/pbarratio;
            pos(2)=(axispos(4)-pos(4))/2+axispos(2);
        else
            pos(2)=axispos(2);
            pos(4)=axispos(4);
            pos(3)=axispos(4)*pbarratio;
            pos(1)=(axispos(3)-pos(3))/2+axispos(1);
        end
    end
end

% convert plot box position to the units used by the axis
% - requires plotting a new axes and removing it
temp=axes('units','pixels','position',pos,...
    'visible','off','parent',get(ax,'parent'));
set(temp,'units',currunit);
pos=get(temp,'position');
delete(temp);

end

