function [list]=available_1dmodels
%AVAILABLE_1DMODELS    Returns available 1D Earth model functions
%
%    Usage:    list=available_1dmodels
%
%    Description:
%     LIST=AVAILABLE_1DMODELS returns a cell array of strings that
%     correspond to functions that provide 1D reference Earth models.  This
%     is just a way of centralizing what is available so I only have to
%     update here (hopefully).
%
%    Notes:
%
%    Examples:
%     % This is mainly for internal usage.  But is kinda
%     % nice to get a quick list of available models:
%     available_1dmodels
%
%    See also: PREM, AK135, IASP91

%     Version History:
%        May  19, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 17:25 GMT

% todo:

% available
list={'ak135' 'iasp91' 'prem'};

end
