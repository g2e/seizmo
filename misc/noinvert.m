function noinvert(fh)
%NOINVERT    Turns off hardcopy black/white inversion
%
%    Usage:    noinvert(fh)
%
%    Description:
%     NOINVERT(FH) turns off hardcopy black & white inversion performed by
%     Matlab.  This is useful for figures with a black background
%     (something I do often here).  Obviously this is for presentations --
%     not for actual printing!
%
%    Notes:
%     - I usually disdain making a function for a one-liner, but this one
%       is used often enough and is annoying enough that I wrote it.
%
%    Examples:
%     % Make and save a SEIZMO figure:
%     fh=plot0(data);
%     fig2print(fh,'landscape');
%     noinvert(fh);
%     print(fh,'-dpdf','tmp.pdf')
%
%    See also: FIG2PRINT

%     Version History:
%        Aug.  4, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check figure handles exist
if(any(~ishandle(fh)) || any(~strcmpi('figure',get(fh,'type'))))
    error('seizmo:noinvert:badHandle',...
        'FH must be a valid figure handle!');
end

% turn off hardcopy color inversion
set(fh,'inverthardcopy','off');

end
