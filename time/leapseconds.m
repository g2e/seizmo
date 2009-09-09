function [leaps]=leapseconds()
%LEAPSECONDS    Returns date and time strings of each leap second
%
%    Usage:    leaps=leapseconds()
%
%    Description: LEAPSECONDS() returns a formatted character array of the
%     date-times of the end of each leap second, which is always the start
%     of the next UTC month.  The last characters of each row gives a '+'
%     or '-' to indicate whether the leap second was positive or negative
%     (all leap seconds have been positive thus far).  For example, the
%     first leap second was positive and occurred on June 30, 1972 at
%     23:59:60 and ended on July 1, 1972 at 00:00:00, so it's string is
%     '01-Jul-1972 00:00:00+'.  A negative leap second would start at the
%     58th second (if we begin each minute with the 0th second), removing
%     one second as there would be no 59th second.
%
%    Notes:
%     - At some point it may be worthwhile to just read in
%       http://maia.usno.navy.mil/ser7/leapsec.dat
%     - Since leapseconds were introduced in 1972, UTC times before that
%       are not properly corrected here to maintain timing near UT1 (GMT).
%       There was actually a different method implemented but that is not
%       a matter for this function.
%     - Q: How to properly handle data around leap seconds from equipment
%          that doesn't handle leap seconds?  
%       A: Leave the last data segment before a leap second alone.  The
%          first and subsequent data segments after a leap second should be
%          shifted to one second earlier until the clock corrected to the
%          actual UTC time (should have to jump back 1 second).  These will
%          then be correct in UTC time.
%     - Q: What if a data segment begins within a leap second on equipment
%          that doesn't handle leap seconds?
%       A: You need to merge that data segment with the prior (without UTC
%          awareness) and then shift all subsequent records to 1 second
%          prior (until the clock locks to the correct UTC time).
%
%    Examples:
%     Ever wondered when those pesky leap seconds were? Just run:
%      leapseconds
%
%    See also: getleapseconds, isleapyear, leapsinday, totalleaps,
%              fixtimes, timediff, utc2tai, tai2utc

%     Version History:
%        Oct. 28, 2008 - initial version
%        Mar. 29, 2009 - added some notes on data handling
%        Apr. 23, 2009 - move usage up
%        June 10, 2009 - minor doc update
%        Sep.  5, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2009 at 19:05 GMT

% todo:

% add dates at the bottom (keep time increasing)
leaps=['01-Jul-1972 00:00:00+'
       '01-Jan-1973 00:00:00+'
       '01-Jan-1974 00:00:00+'
       '01-Jan-1975 00:00:00+'
       '01-Jan-1976 00:00:00+'
       '01-Jan-1977 00:00:00+'
       '01-Jan-1978 00:00:00+'
       '01-Jan-1979 00:00:00+'
       '01-Jan-1980 00:00:00+'
       '01-Jul-1981 00:00:00+'
       '01-Jul-1982 00:00:00+'
       '01-Jul-1983 00:00:00+'
       '01-Jul-1985 00:00:00+'
       '01-Jan-1988 00:00:00+'
       '01-Jan-1990 00:00:00+'
       '01-Jan-1991 00:00:00+'
       '01-Jul-1992 00:00:00+'
       '01-Jul-1993 00:00:00+'
       '01-Jul-1994 00:00:00+'
       '01-Jan-1996 00:00:00+'
       '01-Jul-1997 00:00:00+'
       '01-Jan-1999 00:00:00+'
       '01-Jan-2006 00:00:00+'
       '01-Jan-2009 00:00:00+'];

end

