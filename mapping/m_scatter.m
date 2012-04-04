function [han]=m_scatter(varargin)
%M_SCATTER    Scatter/Bubble plot on an M_MAP map
%
%    Usage: M_SCATTER(LON,LAT)
%           M_SCATTER(LON,LAT,SIZE)
%           M_SCATTER(LON,LAT,SIZE,COLOR)
%           M_SCATTER(...,MARKER)
%           M_SCATTER(...,'filled')
%           M_SCATTER(...,'PropertyName',propertyvalue) 
%           M_SCATTER(axes_handle,...) 
%           H=M_SCATTER(...)
%
%    Description:
%     M_SCATTER provides a wrapper for making a scatter plot on an M_MAP
%     map.  This adjusts the LON & LAT arguments to map coordinates.  Note
%     that the SIZE argument is still in points^2.  All options for SCATTER
%     are available.
%
%    See also: SCATTER

%     Version History:
%        May   9, 2010 - initial version (based on M_PLOT & SCATTER)
%        June 25, 2010 - point clipping to eliminate off-map points
%        Feb. 10, 2011 - minor H1 change
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 22:15 GMT

global MAP_PROJECTION

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

% check for axes, pv pairs
start=1;
if (nargin > 0) && (numel(varargin{1}) == 1) && ...
      ishandle(varargin{1}) && strcmp(get(varargin{1},'type'),'axes')
  start=2;
end
[args,pvpairs]=parseparams(varargin(start:end));
nargs=numel(args);

% are there enough data inputs?
msg=nargchk(2,4,nargs);
if ~isempty(msg);
  disp(msg);
  help m_scatter
  return
end

% convert lon/lat to x/y
% - clips points outside map boundaries
[args{1:2}] = m_ll2xy(args{1:2},'clip','point');

% pass arguments to scatter
if(start==1)
    h=scatter(args{:},pvpairs{:});
else
    h=scatter(varargin{1},args{:},pvpairs{:});
end

if nargout == 1
  han = h;
end

return