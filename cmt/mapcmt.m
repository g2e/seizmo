function [varargout]=mapcmt(mt,varargin)
%MAPCMT    Plots GlobalCMT moment tensors on a M_Map map
%
%    Usage:    mapcmt(mt)
%              mapcmt(mt,'param1',val1,...,'paramN',valN)
%              h=mapcmt(...)
%
%    Description:
%     MAPCMT(MT) creates a map and plots the moment tensors in MT on it.
%     MT must be a scalar struct with a format as that returned by FINDCMTS
%     & FINDCMT.
%
%     MAPCMT(MT,'PARAM1',VAL1,...,'PARAMN',VALN) passes on parameters to
%     PLOTMT to refine how the cmts are plotted.
%
%     H=MAPCMT(...) returns the handles to the moment tensor objects.
%
%    Notes:
%     - CMTs as implimented are incompatible with mmaps because they
%       require shading which breaks mmap patches.  We draw an invisible
%       axis over the map to get around this for now.
%     - The workaround mentioned above still does not allow for saving the
%       figure at vector graphics quality.  Yet another cdata issue.
%     - The overlap of cmt outlines is a matlab bug.  I have turned them
%       off by default here to avoid the annoyance.
%
%    Examples:
%     % 
%
%    See also: FINDCMT, FINDCMTS, PLOTMT, RADPAT, MMAP

%     Version History:
%        May  15, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  15, 2011 at 23:55 GMT

% todo:
% - cmts fail to work with mmaps due to shading issues ... >_<

% check nargin
error(nargchk(1,inf,nargin));

% extract axes handle
[ax,varargin]=axescheck(varargin{:});

% check for map (new if none)
if(isempty(ax))
    % draw one for the lazy
    ax=mmap;
elseif(isscalar(ax) && isreal(ax) ...
        && ishandle(ax) && strcmp('axes',get(ax,'type')))
    axes(ax);
else
    % needs to be a mmap
    error('seizmo:mapcmt:badAxes',...
        'AXHANDLE must be a handle to an appropriate map!');
end

% require map
global MAP_PROJECTION MAP_VAR_LIST
if(isempty(MAP_PROJECTION))
    error('seizmo:mapcmt:badInput',...
        'AXHANDLE must be a handle to an appropriate map!');
end

% hack: lay invisible axis over map with same limits
drawnow;
newax=axes('position',get(ax,'position'),'parent',get(ax,'parent'),...
    'visible','off','color','none','xcolor','w','ycolor','w');
set(newax,'userdata',linkprop([ax newax],{'position' 'xlim' 'ylim' ...
    'parent' 'PlotBoxAspectRatio' 'PlotBoxAspectRatioMode'}));

% set hold on
held=ishold(newax);
hold(newax,'on');

% get lat/lon => x/y
lat=mt.centroidlat;
lon=mt.centroidlon;
[x,y]=m_ll2xy(lon,lat,'clip','point');

% get moment tensor
mt=[mt.mrr mt.mtt mt.mpp mt.mrt mt.mrp mt.mtp]; % Nx6

% get local north
[x1,y1]=m_ll2xy(lon,lat+diff(MAP_VAR_LIST.lats)/1000,'clip','point');
localnorth=90-atan2(y1-y,x1-x)*180/pi;

% let plotmt do the work
h=plotmt(x,y,mt,'roll',localnorth,'r',.1,'lc',[],varargin{:},...
    'parent',newax);
if(nargout); varargout{1}=h; end

% release hold
if(~held); hold(newax,'off'); end

end
