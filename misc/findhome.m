function [home]=findhome()
%FINDHOME    Returns the user's home directory
%
%    Usage:    home=findhome
%
%    Description:
%     HOME=FINDHOME returns the home directory of the user.
%
%    Notes:
%     - Uses environment variables
%
%    Examples:
%     % Go to the home directory:
%     cd(findhome);
%
%    See also: CD, PWD, GETENV

%     Version History:
%        Feb. 16, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2012 at 17:15 GMT

% todo:

% use environment variables
if(ispc); home=getenv('USERPROFILE');
else home=getenv('HOME');
end

end
