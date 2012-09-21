function [lgc]=isgeofss(s)
%ISGEOFSS    True if input is a struct as defined by GEOFSS functions
%
%    Usage:    lgc=isgeofss(s)
%
%    Description:
%     LGC=ISGEOFSS(S) returns TRUE if S is a struct as output by GEOFSS
%     functions and FALSE otherwise.  S may be a nonscalar struct (LGC is
%     still either TRUE or FALSE).  See CHKGEOFSS for more details.
%
%    Notes:
%
%    Examples:
%     % See if GEOFSSXC conforms:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     isgeofss(geofssxc(xcdata,[lat(:) lon(:)],30,[.01 .0125]))
%
%    See also: CHKGEOFSS, CHKFSS, ISFSS

%     Version History:
%        Oct. 10, 2010 - initial version
%        June  4, 2012 - adapted from isgeofkstruct
%        Sep. 12, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
lgc=isempty(chkgeofss(s));

end
