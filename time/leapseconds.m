function [dates,leaps,offsets]=leapseconds(option)
%LEAPSECONDS    Returns leapsecond info (with caching for multiple calls)
%
%    Usage:    [dates,leaps,offsets]=leapseconds()
%              [dates,leaps,offsets]=leapseconds(true|false)
%
%    Description:
%     [DATES,LEAPS,OFFSETS]=LEAPSECONDS() returns leap second info and
%     stores that info in a global variable to speed up subsequent calls.
%     DATES contains the serial dates on which the 'leaps' end, LEAPS
%     contains the corresponding 'leap' (either 1 or -1), and OFFSET gives
%     the corresponding UTC lag from TAI.
%
%     [DATES,LEAPS,OFFSETS]=LEAPSECONDS(OPTION) controls the behavior of
%     the info caching.  OPTION must be TRUE or FALSE (logical).  TRUE will
%     use the cached leapsecond info if there is any otherwise retreiving
%     it via the leapsec.dat included with this mfile.  FALSE will force
%     reading of leapsec.dat to get the leapsecond info. The default is
%     TRUE.
%
%    Notes:
%     - Update the leap second info using LEAPSECONDS_UPDATE
%     - Since leapseconds were introduced in 1972, UTC times before that
%       are not properly accounted for with leap seconds to maintain timing
%       near UT1 (GMT).  There was actually a different method implemented
%       but that is not a matter for this function.  The data for that is
%       found here: http://maia.usno.navy.mil/ser7/tai-utc.dat
%     - Q: How to properly handle data around leap seconds from recording
%          equipment that doesn't handle leap seconds as they occur?  
%       A: Leave the last data segment before a leap second alone.  The
%          first and subsequent data segments after a leap second should be
%          shifted to one second earlier until the clock corrected to the
%          actual UTC time (should have to jump back 1 second).  These will
%          then be correct in UTC time.
%     - Q: What if a data segment begins within a leap second on recording
%          equipment that doesn't handle leap seconds as they occur?
%       A: You need to merge that data segment with the prior (without UTC
%          awareness) and then shift all subsequent records to 1 second
%          prior (until the clock locks to the correct UTC time).  MELD &
%          TIMESHIFT are useful for this task.
%
%    Examples:
%     % Ever wondered when those pesky leap seconds were? Just run:
%     datestr(leapseconds)
%
%     % To see the speed benefit of caching:
%     tic; for i=1:100; [dates,leaps,offsets]=leapseconds(true); end; toc
%     tic; for i=1:100; [dates,leaps,offsets]=leapseconds(false); end; toc
%
%    See also: LEAPSECONDS_UPDATE, UTCOFFSET, LOD, FIXTIMES, TIMEDIFF,
%              UTC2TAI, TAI2UTC, ISLEAPYEAR

%     Version History:
%        Nov. 10, 2008 - initial version
%        Nov. 15, 2008 - update for new name schema, option now a logical
%        Dec. 13, 2008 - fix bug, eliminate recursion
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - minor doc update
%        Aug.  4, 2009 - strictly formatted string passed to datenum
%        Sep.  5, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - renamed from GETLEAPSECONDS to LEAPSECONDS as that
%                        function was moved to a file that can be updated,
%                        added offset from TAI as a 3rd output
%        Feb.  7, 2012 - doc update (merge to meld update)
%        Feb. 13, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% check option
if(nargin==0 || isempty(option))
    option=true;
elseif(~islogical(option))
    error('seizmo:leapseconds:optionBad',...
        'OPTION must be logical!');
end

% global
global SEIZMO

% retrieve from memory or recalculate
if(option && isfield(SEIZMO,'LEAPSECONDS')...
        && isfield(SEIZMO.LEAPSECONDS,'DATES')...
        && isfield(SEIZMO.LEAPSECONDS,'LEAPS')...
        && isfield(SEIZMO.LEAPSECONDS,'OFFSETS')...
        && isnumeric(SEIZMO.LEAPSECONDS.DATES)...
        && isnumeric(SEIZMO.LEAPSECONDS.LEAPS)...
        && isnumeric(SEIZMO.LEAPSECONDS.OFFSETS)...
        && isequal(size(SEIZMO.LEAPSECONDS.DATES),...
        size(SEIZMO.LEAPSECONDS.LEAPS),...
        size(SEIZMO.LEAPSECONDS.OFFSETS)))
    dates=SEIZMO.LEAPSECONDS.DATES;
    leaps=SEIZMO.LEAPSECONDS.LEAPS;
    offsets=SEIZMO.LEAPSECONDS.OFFSETS;
else
    % parse leapsec.dat info
    path=fileparts(mfilename('fullpath'));
    txt=reshape(readtxt([path filesep 'leapsec.dat']),81,[])';
    dates=datenum(txt(:,1:12),'yyyy mmm dd');
    offsets=str2double(cellstr(txt(:,39:42)));
    leaps=[0; diff(offsets)];
    SEIZMO.LEAPSECONDS.DATES=dates;
    SEIZMO.LEAPSECONDS.LEAPS=leaps;
    SEIZMO.LEAPSECONDS.OFFSETS=offsets;
end

end
