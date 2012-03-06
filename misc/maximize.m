function []=maximize(fh)
%MAXIMIZE    Maximize figure in Matlab (uses java, platform independent)
%
%    Usage:    maximize
%              maximize(fh)
%
%    Description:
%     MAXIMIZE maximizes the current figure.
%
%     MAXIMIZE(FH) maximizes the indicated figure.
%
%    Notes:
%     - Does not work with Octave.
%
%    Examples:
%     % Plot and maximize:
%     imagesc(magic(100));
%     maximize;
%
%    See also: SET, GET, FIGURE, HANDLE

%     Version History:
%        Mar.  6, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2012 at 11:15 GMT

% todo:

% use current figure if none indicated
if(nargin==0); fh=gcf; end

% assure figure is drawn
drawnow;

% maximize through java
obj=get(handle(fh),'JavaFrame'); 
obj.setMaximized(true);

end
