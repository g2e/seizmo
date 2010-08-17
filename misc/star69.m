function [who]=star69(n)
%STAR69    Returns who called the current function
%
%    Usage:  caller=star69
%            caller=star69(gen)
%
%    Description:
%     CALLER=STAR69 returns the caller of the function that has called this
%     function.  That is, it returns the grandparent of STAR69.  So if FUN1
%     calls FUN2 and FUN2 uses STAR69, STAR69 will return FUN1.  If there
%     is no calling function (as in it was directly called from the command
%     line) then CALLER is set to ''.
%
%     CALLER=STAR69(GEN) specifies which caller generation to return.  The
%     default is 1 (parent).  A value of 0 will return the current function
%     and a value of 2 will return the grandparent (if it exists).
%
%    Notes:
%
%    Examples:
%     % grandparent (2 generations up) of current function
%     star69(2)
%
%    See also: DBSTACK, MFILENAME

%     Version History:
%        Aug. 14, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 14, 2010 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default/check generation
if(nargin<1 || isempty(n)); n=1; end
if(~isreal(n) || ~isscalar(n) || n~=fix(n))
    error('seizmo:star69:badInput',...
        'GEN must be a scalar integer!');
end
n=n+2;

% get calling stack
% 1 is star69
% 2 is star69's caller
% 3 is star69's caller's caller (default)
callers=dbstack;

% get caller
if(n<1 || n>numel(callers))
    who='';
else
    who=callers(n).name;
end

end
