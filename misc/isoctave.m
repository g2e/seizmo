function [tf]=isoctave()
%ISOCTAVE    Returns TRUE if the application is Octave
%
%    Usage:    tf=isoctave
%
%    Description:
%     ISOCTAVE returns TRUE if this is run in GNU Octave or returns FALSE
%     if this is run in Matlab.
%
%    Notes:
%
%    Examples:
%     % Useful for saving:
%     if(isoctave)
%         save -7 blah.mat myvar
%     else % matlab
%         save blah.mat myvar
%     end
%
%    See also: GETAPPLICATION, VER, VERSION

%     Version History:
%        Mar.  1, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2012 at 19:55 GMT

% todo:

% test for octave interpreter
if(exist('OCTAVE_VERSION','builtin')==5)
    tf=true;
else
    tf=false;
end

end
