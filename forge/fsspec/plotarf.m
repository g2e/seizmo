function [varargout]=plotarf(s,varargin)
%PLOTARF    Plots ARF spectra with source slownesses marked
%
%    Usage:    plotarf(s)
%              plotarf(s,...)
%              ax=plotarf(...)
%
%    Description:
%     PLOTARF(S) plots the ARF (a frequency-slowness power spectra) in
%     struct S.  See ARF for details on S.  This just calls PLOTFSS and
%     then marks the source slownesses afterwards with stars.
%
%     PLOTARF(S,...) passes additional options to PLOTFSS.
%
%     AX=PLOTARF(...) returns the handle to the axis plotted in.
%
%    Notes:
%     - The handle to the stars is tagged as 'ARF_source'.
%
%    Examples:
%     % Get a multi-plane wave response at 0.03Hz for a random array:
%     s=arf(randlatlon(20)/45,50,201,[0 270 45],30,0.03);
%     plotarf(fssavg(s));
%
%    See also: ARF, PLOTFSS, FSSAVG, FSSSUB

%     Version History:
%        June 12, 2012 - initial version
%        Sep. 15, 2012 - adapted from plotgeoarf
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 15, 2012 at 15:05 GMT

% todo:

error(nargchk(1,inf,nargin));
ax=plotfss(s,varargin{:});
userdata=get(ax,'userdata');
hold(ax,'on');
x=s.source.slow.*sind(s.source.baz);
y=s.source.slow.*cosd(s.source.baz);
plot(ax,x,y,'p',...
    'markeredgecolor',userdata.bgcolor,...
    'tag','ARF_source');
hold(ax,'off');
if(nargout); varargout={ax}; end

end
