function [varargout]=plotgeoarf(s,varargin)
%PLOTGEOARF    Plots geoARF spectra with source locations marked
%
%    Usage:    plotgeoarf(s)
%              plotgeoarf(s,...)
%
%    Description:
%     PLOTGEOARF(S) plots the geoARF (a frequency-slowness-position power
%     spectra) in struct S.  See GEOARF for details on S.  This just calls
%     PLOTGEOFSS and then plots the source locations afterwards.
%
%     PLOTGEOARF(S,...) passes additional options to PLOTGEOFSS.
%
%    Notes:
%
%    Examples:
%     % Global ARF for a random array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,[0 0],25,1/500);
%     plotgeoarf(s);
%
%     % Multi-source ARF for a global array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,randlatlon(5),25,1/500);
%     plotgeoarf(geofssavg(s));
%
%    See also: GEOARF, PLOTGEOFSS, GEOFSSAVG, GEOFSSSUB

%     Version History:
%        June 12, 2012 - initial version
%        Sep. 29, 2012 - handle extended lalo
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 15:05 GMT

% todo:

error(nargchk(1,inf,nargin));
ax=plotgeofss(s,varargin{:});
hold(ax,'on');
% have to do this afterwards b/c of matlab warnings/issues
mmap('ev',s.source.latlon(:,1:2),'parent',ax,varargin{3:end});
hold(ax,'off');
if(nargout); varargout={ax}; end

end
