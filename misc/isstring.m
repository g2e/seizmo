function [lgc]=isstring(str)
%ISSTRING    True for a string (row vector) of characters
%
%    Usage:    lgc=isstring(str)
%
%    Description:
%     LGC=ISSTRING(STR) returns TRUE if STR is a string (ie row vector) of
%     characters.  This means ISCHAR(STR) must be TRUE and SIZE(STR,1)==1.
%
%    Notes:
%
%    Examples:
%     % A 2x2 character array will return FALSE:
%     isstring(repmat('a',2,2))
%
%    See also: ISCHAR

%     Version History:
%        Sep. 13, 2010 - initial version (added docs)
%        Oct. 11, 2010 - allow empty string
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2010 at 11:00 GMT

lgc=ischar(str) && ndims(str==2) ...
    && (size(str,1)==1 || isequal(size(str),[0 0]));

end