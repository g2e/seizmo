function [list]=available_3dmodels
%AVAILABLE_3DMODELS    Returns available 3D Earth model functions
%
%    Usage:    list=available_3dmodels
%
%    Description: LIST=AVAILABLE_3DMODELS returns a cell array of strings
%     that correspond to functions that provide 3D reference Earth models.
%     This is just a way of centralizing what is available so I only have
%     to update here (hopefully).
%
%    Notes:
%
%    Examples:
%     This is mainly for internal usage.  But is kinda nice to get a quick
%     list of available models:
%      available_3dmodels
%
%    See also: GET_SCRIPPS_POINT, HMSL06S, HMSL06P, SB4L18

%     Version History:
%        May  19, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  19, 2010 at 17:25 GMT

% todo:

% available
list={'HMSL06S' 'HMSL06P' 'SB4L18'};

end
