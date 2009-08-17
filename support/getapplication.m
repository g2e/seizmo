function [application,version]=getapplication()
%GETAPPLICATION    Returns application running this script and its version
%
%    Usage:    [application,version]=getapplication()
%
%    Description: [APPLICATION,VERSION]=GETAPPLICATION() will determine and
%     return the name and version of the application running this script
%     (obviously only if the application can run this script in the first
%     place).  Both APPLICATION and VERSION are strings.
%
%    Notes:
%     - returns 'UNKNOWN' if it cannot figure out the application
%
%    Examples:
%     Matlab and Octave still behave quite differently for a number of
%     different functions so it is best in some cases to use different
%     function calls depending on which we are running:
%      [app,ver]=getapplication;
%      if(strcmp(app,'MATLAB'))
%        % do something via matlab routines
%      else
%        % do something via octave routines
%      end
%
%    See also:

%     Version History:
%        Nov. 13, 2008 - initial version
%        Mar.  3, 2009 - minor doc cleaning
%        Apr. 23, 2009 - move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% checking for Matlab will throw an error in Octave
try
    % first check if we are in Matlab
    a=ver('matlab');
    
    % we are in Matlab
    application=a.Name;
    version=a.Version;
    return;
catch
    % check if we are in Octave
    if(exist('OCTAVE_VERSION','builtin')==5)
        application='OCTAVE';
        version=OCTAVE_VERSION;
        return;
    % ok I have no clue what is running
    else
        application='UNKNOWN';
        version='UNKNOWN';
        return;
    end
end

end
