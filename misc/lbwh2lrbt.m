function [lrbt]=lbwh2lrbt(lbwh)
%LBWH2LRBT    Convert left-bottom-width-height to left-right-bottom-top
%
%    Usage:    lrbt=lbwh2lrbt(lbwh)
%
%    Description:
%     LRBT=LBWH2LRBT(LBWH) converts left/bottom/width/height position
%     coordinates to left/right/bottom/top.  This is useful for manual axes
%     position setting and getting margins.
%
%    Notes:
%
%    Examples:
%     % Find the margins of a axis:
%     lbwh2lrtb(get(gca,'position'))
%
%    See also: LRBT2LBWH

%     Version History:
%        Aug.  4, 2010 - initial version
%        Apr. 13, 2011 - fix h1 line
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check input
if(~isreal(lbwh) || ndims(lbwh)~=2 || size(lbwh,2)~=4)
    error('seizmo:lbwh2lrbt:badInput',...
        'LBWH must be an Nx4 array of [Left Bottom Width Height]');
end

% convert
lrbt=[lbwh(1) lbwh(1)+lbwh(3) lbwh(2) lbwh(2)+lbwh(4)];

end
