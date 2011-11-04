function [caldates]=doy2cal(doydates)
%DOY2CAL    Convert Year & Day-of-Year to Year & Month & Day-of-Month
%
%    Usage:    caldates=doy2cal(doydates)
%
%    Description:
%     CALDATES=DOY2CAL(DOYDATES) returns a Nx3 numeric array of calendar
%     dates CALDATES corresponding to DOYDATES, a Nx2 numeric array with
%     the first column having the year and the second the day-of-year.  The
%     output array has the first column as the year, the second as the
%     month number (1-12) and the last column with the day-of-month.
%     Handles leapdays.  Gregorian calendar only.
%
%    Notes:
%
%    Examples:
%     % Get the calendar date of today:
%     caldate=doy2cal([0 now])
%
%    See also: CAL2DOY, ISLEAPYEAR, FIXDATES, FIXTIMES, TIMEDIFF,
%              GREGORIAN2MODSERIAL, GREGORIAN2SERIAL, SERIAL2GREGORIAN,
%              MODSERIAL2GREGORIAN

%     Version History:
%        Oct. 31, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        June 10, 2010 - handle empty
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 23:45 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require Nx2 numeric
sz=size(doydates);
if(~isnumeric(doydates) || sz(2)~=2)
    error('seizmo:doy2cal:badInput',...
        'DATES array must have Nx2 elements: [years dayofyears]');
elseif(sz(1)==0)
    caldates=zeros(0,3,class(doydates));
    return;
end

% fix dates
doydates=fixdates(doydates);

% permute to make things easier
ndates=prod(sz)/2;
doydates=permute(doydates,[2 1 3:numel(sz)]);
caldates=zeros([3 sz(1) sz(3:end)],class(doydates));

% year is good
caldates(1,:)=doydates(1,:);

% which are leap years
leap=isleapyear(doydates(1,:)).';

% doy of the start of months
ndays=cumsum([0,31,28,31,30,31,30,31,31,30,31,30])+1;
ndays=ndays(ones(ndates,1),:);
ndays(:,3:end)=ndays(:,3:end)+leap(:,ones(1,10));

% calendar month and day
caldates(2,:)=sum(doydates(2*ones(12,1),:).'>=ndays,2);
caldates(3,:)=doydates(2,:)-ndays((1:ndates)+(caldates(2,:)-1)*ndates)+1;

% permute back
caldates=permute(caldates,[2 1 3:numel(sz)]);

end
