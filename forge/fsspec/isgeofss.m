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
%     isgeofss(geofssxc(xcdata,latlon,slow,frng))
%
%    See also: CHKGEOFSS, CHKFSS, ISFSS

%     Version History:
%        Oct. 10, 2010 - initial version
%        June  4, 2012 - adapted from isgeofkstruct
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  4, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
lgc=isempty(chkgeofss(s));

end
