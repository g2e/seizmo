function [leaps]=leapseconds()
%LEAPSECONDS    Returns date and time strings of each leap second
%
%    Description: LEAPSECONDS() returns a formatted character array of the
%     date-times of the end of each leap second, which is always the start
%     of the next UTC month.  The last characters of each row gives a '+'
%     or '-' to indicate whether the leap second was positive or negative
%     (all leap seconds have been positive thus far).  For example, the
%     first leap second was positive and occurred on June 30, 1972 at
%     23:59:60 and ended on July 1, 1972 at 00:00:00, so it's string is
%     '01-Jul-1972 00:00:00+'.  A negative leap second would start at the
%     58th second, removing one second as there would be no 59th second.
%
%    Notes:
%     - At some point it may be worthwhile to just read in
%       http://maia.usno.navy.mil/ser7/leapsec.dat
%     - Since leapseconds were introduced in 1972, UTC times before that
%       are not properly corrected here to maintain timing near UT1 (GMT).
%       There was actually a different method implemented but that is not
%       a matter for this function.
%
%    Tested on: Matlab r2007b
%
%    Usage:    leaps=leapseconds()
%
%    Examples:
%     Ever wondered when those pesky leap seconds were? Just run:
%      leapseconds
%
%    See also: julday, calday, isleapyear

%     Version History:
%        Oct. 28, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 28, 2008 at 15:00 GMT

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

