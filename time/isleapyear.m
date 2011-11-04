function [isleap]=isleapyear(year)
%ISLEAPYEAR    True if year is a leap year
%
%    Usage:    leapyears=isleapyear(years)
%
%    Description:
%     LEAPYEARS=ISLEAPYEAR(YEARS) returns a logical array equal in size to
%     YEARS with values set to true for the corresponding elements in YEARS
%     that are leap years.
%
%    Notes:
%     - Gregorian calendar only!
%     - Only valid for 1+ AD unless corrections are made to BC years
%       (1 BC == 0, 2 BC == -1, ...)
%
%    Examples:
%     % Takes into account the century rule and the exception so that
%     isleapyear([1900 1904 2000 2004]) returns [0 1 1 1].
%
%    See also: FIXDATES, FIXTIMES, TIMEDIFF, CAL2DOY, DOY2CAL,
%              GREGORIAN2MODSERIAL, GREGORIAN2SERIAL, SERIAL2GREGORIAN,
%              MODSERIAL2GREGORIAN

%     Version History:
%        Oct. 28, 2008 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 19:10 GMT

% todo:

% require numeric years
if(~isnumeric(year))
    error('seizmo:isleapyear:badYear','YEARS must be numeric!');
end

% force years to be integer
year=floor(year);

% every 4 years unless it falls on a century except every 400 years
isleap=(~mod(year,4) & (mod(year,100) | ~mod(year,400)));

end
