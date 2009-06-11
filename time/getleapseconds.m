function [dates,leaps]=getleapseconds(option)
%GETLEAPSECONDS    Returns leapsecond info with caching for multiple calls
%
%    Usage:    [dates,leaps]=getleapseconds()
%              [dates,leaps]=getleapseconds(true|false)
%
%    Description:  [DATES,LEAPS]=GETLEAPSECONDS() returns leap second info
%     and stores that info in a global variable to speed up subsequent
%     calls.  DATES contains the serial dates on which the 'leaps' end
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
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - minor doc update, add testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2009 at 00:00 GMT

% todo:

% check nargin
msg=nargchk(0,1,nargin);
if(~isempty(msg)); error(msg); end;

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
    leaps=2*strcmp('+',cellstr(leapstr(:,end)))-1; % +/-1
    [dates,idx]=sort(dates);
    leaps=leaps(idx);
    SEIZMO.GETLEAPSECONDS.DATES=dates;
    SEIZMO.GETLEAPSECONDS.LEAPS=leaps;
end

end
