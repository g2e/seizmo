function [lgc]=isfss(s)
%ISFSS    True if input is a struct as defined by FSS functions
%
%    Usage:    lgc=isfss(s)
%
%    Description:
%     LGC=ISFSS(S) returns TRUE if S is a struct as output by FSS
%     functions and FALSE otherwise.  S may be a nonscalar struct (LGC is
%     still either TRUE or FALSE).  See CHKFSS for more details.
%
%    Notes:
%
%    Examples:
%     % See if FSSXC conforms:
%     isfss(fssxc(xcdata,50,101,[.01 .0125]))
%
%    See also: CHKFSS, CHKGEOFSS, ISGEOFSS

%     Version History:
%        Oct. 10, 2010 - initial version
%        June  4, 2012 - adapted from isgeofkstruct
%        Sep. 12, 2012 - adapted from isgeofss
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
lgc=isempty(chkfss(s));

end
