function [varargout]=mapcmts(cmts,varargin)
%MAPCMTS    Plots GlobalCMT moment tensors on a M_Map map
%
%    Usage:    mapcmts(cmts)
%              mapcmts(cmts,'param1',val1,...,'paramN',valN)
%              h=mapcmts(...)
%
%    Description:
%     MAPCMTS(CMTS) creates a map with the moment tensors in CMTS on it.
%     CMTS must be a scalar struct with a format as that returned by
%     FINDCMTS & FINDCMT.  Note that the beachballs are rotated to point
%     northward instead of just to 'up' on the plot.
%
%     MAPCMTS(CMTS,'PARAM1',VAL1,...,'PARAMN',VALN) passes parameters on to
%     PLOTMT to refine how the moment tensors are plotted.
%
%     H=MAPCMTS(...) returns the handles to the moment tensor objects.
%
%    Notes:
%     - Matlab has trouble keeping the order of many patches so the plot
%       may look bad (the moment tensor patches will be hidden behind the
%       continents and the outlines will be on top of all patches) if there
%       are many moment tensors (around 100 or more).  The moment tensors
%       will be drawn correctly in a pdf export though (see EXPORT_FIG).
%
%    Examples:
%     % Map the first 100 moment tensors in the GlobalCMT catalog:
%     mapcmts(findcmt('n',100));
%
%     % Plot the biggest 100 moment tensors in the GlobalCMT catalog:
%     mapcmts(findcmt('n',100,'magnitude',10));
%
%    See also: FINDCMT, FINDCMTS, PLOTMT, RADPAT, MMAP, EXPORT_FIG

%     Version History:
%        May  15, 2011 - initial version
%        May  31, 2011 - work with contour version of plotmt
%        Jan. 11, 2012 - rename from mapcmt to mapcmts, edit note
%        Feb.  7, 2012 - minor doc update
%        Mar. 23, 2013 - doc update, minor code fixes
%        Jan. 27, 2014 - use axparse instead of axescheck for octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% extract axes handle
[ax,varargin]=axparse(varargin{:});

% check cmts is a scalar struct with required fields
valid={'centroidlat' 'centroidlon' 'mrr' 'mtt' 'mpp' 'mrt' 'mrp' 'mtp'};
if(~isstruct(cmts) || ~isscalar(cmts))
    error('seizmo:mapcmts:badInput',...
        'CMTS must be a scalar struct as from FINDCMTS!');
elseif(~all(ismember(valid,fieldnames(cmts))))
    error('seizmo:mapcmts:badInput',...
        'CMTS is must have the following fields:');
end

% check for map (new if none)
if(isempty(ax))
    % draw one for the lazy
    ax=mmap;
elseif(isscalar(ax) && isreal(ax) ...
        && ishandle(ax) && strcmp('axes',get(ax,'type')))
    axes(ax);
else
    % needs to be a mmap
    error('seizmo:mapcmts:badAxes',...
        'AXHANDLE must be a handle to an appropriate map!');
end

% require map
global MAP_PROJECTION MAP_VAR_LIST
if(isempty(MAP_PROJECTION))
    error('seizmo:mapcmts:badInput',...
        'AXHANDLE must be a handle to an appropriate map!');
end

% set hold on
held=ishold(ax);
hold(ax,'on');

% get lat/lon => x/y
lat=cmts.centroidlat;
lon=cmts.centroidlon;
[x,y]=m_ll2xy(lon,lat,'clip','point');

% get moment tensor
cmts=[cmts.mrr cmts.mtt cmts.mpp cmts.mrt cmts.mrp cmts.mtp]; % Nx6

% get local north (not atomic!)
[x1,y1]=m_ll2xy(lon,lat+diff(MAP_VAR_LIST.lats)/1000,'clip','point');
localnorth=90-atan2(y1-y,x1-x)*180/pi;

% let plotmt do the work
h=plotmt(x,y,cmts,'roll',localnorth,'r',.1,varargin{:},'parent',ax);

% output if desired
if(nargout); varargout{1}=h; end

% release hold
if(~held); hold(ax,'off'); end

end
