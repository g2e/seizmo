function [lbwh]=lrbt2lbwh(lrbt)
%LRBT2LBWH    Convert left-right-bottom-top to left-bottom-width-height
%
%    Usage:    lbwh=lrbt2lbwh(lrbt)
%
%    Description:
%     LBWH=LRBT2LBWH(LRBT) convert left/right/bottom/top position
%     coordinates to left/bottom/width/height.  This is useful when doing
%     manual specifying of axes positions.  Note that the width and height
%     may be negative if the coordinate system does not increase when going
%     up and/or right.
%
%    Notes:
%
%    Examples:
%     % Change an axis to have .13 padding on each side:
%     set(gca,'position',lrbt2lbwh([.13 .13 .13 .13]));
%
%    See also: LBWH2LRBT

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
if(~isreal(lrbt) || ndims(lrbt)~=2 || size(lrbt,2)~=4)
    error('seizmo:lrbt2lbwh:badInput',...
        'LRBT must be an Nx4 array of [Left Right Bottom Top]');
end

% convert
lbwh=[lrbt(1) lrbt(3) lrbt(2)-lrbt(1) lrbt(4)-lrbt(3)];

end
