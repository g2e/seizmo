function []=raisefancy(ax)
%RAISEFANCY    Raise fancy brackets to front of plot
%
%    Usage:    raisefancy
%              raisefancy(ax)
%
%    Description:
%     RAISEFANCY raises the fancy axes of the current axes if they exist.
%
%     RAISEFANCY(AX) raises the fancy axes of the axes specified by axes
%     handle AX.
%
%    Notes:
%
%    Examples:
%     % Make a map with fancy axes, add some features and raise the fancy
%     % axes above those features:
%     ax=mmap('po',{'lat',[-10 20],'lon',[-10 20]},'go',{'box','fancy'});
%     mapfeature(ax,{'boundaries' 'isochrons' 'fz' 'seamounts' 'lips' ...
%                    'wars' 'cvl' 'jos' 'cc' 'impacts' 'volcanoes' ...
%                    'hotspots'});
%     raisefancy(ax);
%
%    See also: MMAP, MAPFEATURE, MOVEKIDS, FINDOBJ

%     Version History:
%        Feb.  9, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2011 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default is current axes if none specified
if(nargin<1 || isempty(ax)); ax=gca; end

% raise fancy axes if there
h=findobj(ax,'-regexp','tag','m_grid_fancybox*');
if(~isempty(h)); movekids(h,'front');

end
