function [doydates]=cal2doy(caldates)
%CAL2DOY    Convert Year & Month & Day-of-Month to Year & Day-of-Year
%
%    Usage:    doydates=cal2doy(caldates)
%
%    Description:
%     DOYDATES=CAL2DOY(CALDATES) returns a Nx2 numeric array of day-of-year
%     dates corresponding to CALDATES, a Nx3 numeric array with the 1st
%     column giving the year, the 2nd column lists the month number (1-12)
%     and the 3rd column has the day-of-month.  The output array has the
%     1st column as the year and the 2nd column lists the day-of-year.
%     Handles leapdays.  Gregorian calendar only.
%
%    Notes:
%
%    Examples:
%     % Get the day of year of now:
%     doydate=cal2doy([0 1 now])
%
%    See also: DOY2CAL, ISLEAPYEAR, FIXDATES, FIXTIMES, TIMEDIFF,
%              GREGORIAN2MODSERIAL, GREGORIAN2SERIAL, SERIAL2GREGORIAN,
%              MODSERIAL2GREGORIAN

%     Version History:
%        Oct. 31, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require Nx3 numeric
sz=size(caldates);
if(~isnumeric(caldates) || sz(2)~=3)
    error('seizmo:cal2doy:badInput',...
        'DATES array must have Nx3 elements: [years months dayofmonths]');
end

% fix dates
caldates=fixdates(caldates);

% permute to make things easier
ndates=prod(sz)/3;
caldates=permute(caldates,[2 1 3:numel(sz)]);
doydates=zeros([2 sz(1) sz(3:end)]);

% year is good
doydates(1,:)=caldates(1,:);

% which are leap years
leap=isleapyear(caldates(1,:)).';

% doy of the start of months minus 1
ndays=cumsum([0,31,28,31,30,31,30,31,31,30,31,30]);
ndays=ndays(ones(ndates,1),:);
ndays(:,3:end)=ndays(:,3:end)+leap(:,ones(1,10));

% day of year
doydates(2,:)=caldates(3,:)+ndays((1:ndates)+(caldates(2,:)-1)*ndates);

% permute back
doydates=permute(doydates,[2 1 3:numel(sz)]);

end
