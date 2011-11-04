function [dates]=fixdates(dates)
%FIXDATES    Assures calendar or day-of-year dates are in proper range
%
%    Usage:    dates=fixdates(dates)
%
%    Description:
%     FIXDATES(DATES) returns equivalent dates to the dates in DATES such
%     that they follow typical Gregorian calendar conventions (months are 1
%     to 12, dayofmonths are 1 to the last day of the month, and dayofyears
%     are from 1 to 365 or 366 if it is a leap year).  DATES must be a Nx2
%     array of [years dayofyears] or a Nx3 array of
%     [years months dayofmonths].  Inputs should only be integers.
%     Gregorian calendar only!
%
%    Notes:
%
%    Examples:
%     % Nifty way to say the last day of February:
%     fixdates([someyear 3 0])
%     % will return [someyear 2 28] or [someyear 2 29].
%
%     % Five hundred days after the start of some year:
%     fixdates([someyear 500])
%
%    See also: FIXTIMES, DOY2CAL, CAL2DOY, TIMEDIFF, ISLEAPYEAR,
%              GREGORIAN2MODSERIAL, GREGORIAN2SERIAL, SERIAL2GREGORIAN,
%              MODSERIAL2GREGORIAN

%     Version History:
%        Nov.  1, 2008 - initial version
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

% require numeric
if(~isnumeric(dates))
    error('seizmo:fixdates:badInput','DATES must be a numeric array!');
end

% force integers
dates=floor(dates);

% check input
sz=size(dates);
dates=permute(dates,[2 1 3:numel(sz)]);
switch sz(2)
    % yr, doy
    case 2
        % deal with low/high doy values
        low=dates(2,:)<1;
        while(any(low))
            dates(1,low)=dates(1,low)-1;
            dates(2,low)=dates(2,low)+365+isleapyear(dates(1,low));
            low(low)=dates(2,low)<1;
        end
        leap=isleapyear(dates(1,:));
        high=dates(2,:)>(365+leap);
        while(any(high))
            dates(1,high)=dates(1,high)+1;
            dates(2,high)=dates(2,high)-leap(high)-365;
            leap(high)=isleapyear(dates(1,high));
            high(high)=dates(2,high)>(365+leap(high));
        end
    % yr, mon, dom
    case 3
        % deal with low/high months
        dates(1,:)=dates(1,:)+floor((dates(2,:)-1)/12);
        dates(2,:)=mod(dates(2,:)-1,12)+1;
        
        % deal with low/high days
        dim=[31 28 31 30 31 30 31 31 30 31 30 31];
        leap=isleapyear(dates(1,:));
        low=find(dates(3,:)<1);
        while(~isempty(low))
            dates(2,low)=dates(2,low)-1;
            rollover=dates(2,low)<1;
            dates(1,low(rollover))=dates(1,low(rollover))-1;
            leap(low(rollover))=isleapyear(dates(1,low(rollover)));
            dates(2,low(rollover))=12;
            dates(3,low)=dates(3,low)...
                +(dim(dates(2,low))+(leap(low) & dates(2,low)==2));
            low=low(dates(3,low)<1);
        end
        high=find(dates(3,:)>(dim(dates(2,:))+(leap & dates(2,:)==2)));
        while(~isempty(high))
            dates(3,high)=dates(3,high)...
                -(dim(dates(2,high))+(leap(high) & dates(2,high)==2));
            dates(2,high)=dates(2,high)+1;
            rollover=dates(2,high)>12;
            dates(1,high(rollover))=dates(1,high(rollover))+1;
            leap(high(rollover))=isleapyear(dates(1,high(rollover)));
            dates(2,high(rollover))=1;
            high=high(...
                dates(3,high)>(dim(dates(2,high))...
                +(leap(high) & dates(2,high)==2)));
        end
    % unknown
    otherwise
        error('seizmo:fixdates:badInput',...
            'DATES must be Nx2 or Nx3!');
end
dates=permute(dates,[2 1 3:numel(sz)]);

end
