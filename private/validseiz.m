function [valid]=validseiz(filetype)
%VALIDSEIZ    Returns valid SAClab datafile versions
%
%    Description: VALIDSEIZ(FILETYPE) returns a vector of version numbers
%     of the specified filetype FILETYPE that SAClab can work with.  If the
%     filetype is not supported, VALIDSEIZ will return an empty array.
%
%    Notes:
%     - currently only supports filetype 'SAClab Binary'
%     - version 6 corresponds to SAC's v6 binary format
%     - versions 101,200,201 are modifications of SAC's v6.
%
%    Tested on: Matlab r2007b
%    
%    Usage:    valid_versions=validseiz(filetype)
%
%    Examples:
%     How many different SAClab binary file versions are supported:
%      length(validseiz('SAClab Binary'))
%
%    See also:  seizchk, isseiz, seizdef, getversion

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        Apr. 18, 2008 - drop a couple deprecated versions
%        June 12, 2008 - doc update, added history
%        Sep. 14, 2008 - minor doc update
%        Oct. 17, 2008 - filetype argument, history fix
%        Oct. 27, 2008 - minor doc update
%        Nov. 13, 2008 - renamed from VVSEIS to VALIDSEIZ
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2008 at 04:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin))

% SAClab binary version
if(strcmpi(filetype,'SAClab Binary'))
    valid=[6 101 200 201];
end

end
