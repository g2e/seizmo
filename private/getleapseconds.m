function [dates,leaps]=getleapseconds(option)
%GETLEAPSECONDS    Returns leapsecond info with caching for multiple calls
%
%    Description:  [DATES,LEAPS]=GETLEAPSECONDS() returns leap second info
%     and stores that info in a global variable to speed up subsequent
%     calls.  DATES contains the serial dates on which the 'leaps' occur
%     and LEAPS contains the offset (either 1 or -1).
%
%     [DATES,LEAPS]=GETLEAPSECONDS(OPTION) controls the behavior of the
%     info caching.  OPTION must be TRUE or FALSE (logical).  TRUE will use
%     the cached leapsecond info if there is any otherwise retreiving it
%     via LEAPSECONDS.  FALSE will always use LEAPSECONDS to get the latest
%     leapsecond info.  The default is TRUE.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    [dates,leaps]=getleapseconds()
%              [dates,leaps]=getleapseconds(true|false)
%
%    Examples:
%     To see the speed benefit of caching:
%      tic; for i=1:100; [dates,leaps]=getleapseconds(true); end; toc
%      tic; for i=1:100; [dates,leaps]=getleapseconds(false); end; toc
%
%    See also: leapseconds, totalleaps, leapsinday

%     Version History:
%        Nov. 10, 2008 - initial version
%        Nov. 15, 2008 - update for new name schema, option now a logical
%        Dec. 13, 2008 - fix bug, eliminate recursion
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 13, 2008 at 07:50 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% check option
if(nargin==0 || isempty(option))
    option=true;
elseif(~islogical(option))
    error('seizmo:getleapseconds:optionBad',...
        'OPTION must be logical!');
end

% global
global SEIZMO

% retrieve from memory or recalculate
if(option && isfield(SEIZMO,'GETLEAPSECONDS')...
        && isfield(SEIZMO.GETLEAPSECONDS,'DATES')...
        && isfield(SEIZMO.GETLEAPSECONDS,'LEAPS')...
        && isnumeric(SEIZMO.GETLEAPSECONDS.DATES)...
        && isnumeric(SEIZMO.GETLEAPSECONDS.LEAPS)...
        && isequal(size(SEIZMO.GETLEAPSECONDS.LEAPS),...
        size(SEIZMO.GETLEAPSECONDS.DATES)))
    dates=SEIZMO.GETLEAPSECONDS.DATES;
    leaps=SEIZMO.GETLEAPSECONDS.LEAPS;
else
    leapstr=leapseconds;
    dates=datenum(leapstr);
    leaps=2*strcmp('+',cellstr(leapstr(:,end)))-1;
    [dates,idx]=sort(dates);
    leaps=leaps(idx);
    SEIZMO.GETLEAPSECONDS.DATES=dates;
    SEIZMO.GETLEAPSECONDS.LEAPS=leaps;
end

end
