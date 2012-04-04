function [list]=available_3dmodels
%AVAILABLE_3DMODELS    Returns available 3D Earth model functions
%
%    Usage:    list=available_3dmodels
%
%    Description:
%     LIST=AVAILABLE_3DMODELS returns a cell array of strings that
%     correspond to functions that provide 3D reference Earth models.  This
%     is just a way of centralizing what is available so I only have to
%     update here (hopefully).
%
%    Notes:
%
%    Examples:
%     % This is mainly for internal usage.  But is kinda
%     % nice to get a quick list of available models:
%     available_3dmodels
%
%    See also: GET_SCRIPPS_VALUE, DZ04, HMSL06S, HMSL06P, PRI05, S20RTS,
%              SAW24B16, SB4L18, TX2006, TX2007, MIT08

%     Version History:
%        May  19, 2010 - initial version
%        July 25, 2010 - added several models
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 17:25 GMT

% todo:

% available
list={'DZ04' 'HMSL06S' 'HMSL06P' 'MIT-P08' 'PRI-P05' 'PRI-S05' 'S20RTS' ...
      'SAW24B16' 'SB4L18' 'TX2006' 'TX2007'};

end
