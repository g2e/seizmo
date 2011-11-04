function [modserial]=gregorian2modserial(times)
%GREGORIAN2MODSERIAL    Convert Gregorian dates to modified serial dates
%
%    Usage:    modserialdates=gregorian2modserial(gregoriandates)
%
%    Description:
%     [MODSERIAL]=GREGORIAN2MODSERIAL(DATES) returns the equivalent
%     modified serial dates for the Gregorian dates stored in DATES.  A
%     modified serial date of 1 day and 0 seconds corresponds to 
%     January 1, 0000 at 00:00:00 and is in the format [days seconds].
%     DATES must be a Nx2, Nx3, Nx5 or Nx6 numeric array of either 
%     [year dayofyear], [year month dayofmonth],
%     [year dayofyear hour minute seconds] or
%     [year month dayofmonth hour minute seconds].
%
%    Notes:
%     - Do not use partial years or months as they will not be handled
%       correctly!
%     - Does not account for UTC leap seconds (ie assumes everyday has
%       86400 seconds)!
%     - Basically like Matlab's DATENUM except that it also accepts day of
%       year style dates, rolls-over negative calendar months, does not
%       handle string input, and outputs day and secondofday rather than
%       just days (for better precision).
%
%    Examples:
%     % Why use modserial over serial? Better precision:
%     modserial=gregorian2modserial([0 now])
%     serial=gregorian2serial([0 now])
%
%    See also: SERIAL2GREGORIAN, GREGORIAN2SERIAL, MODSERIAL2GREGORIAN,
%              FIXTIMES, FIXDATES, TIMEDIFF, ISLEAPYEAR, CAL2DOY, DOY2CAL

%     Version History:
%        Nov. 11, 2008 - initial version
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

% check times
sz=size(times);
if(~isnumeric(times) || ~any(sz(2)==[2 3 5 6]))
    error('seizmo:gregorian2modserial:badInput',...
        'DATES must be a numeric Nx2, Nx3, Nx5 or Nx6 array!');
end

% work by date type
switch sz(2)
    case 2
        % year, day of year
        days=times(:,1,:)*365+floor(times(:,2,:));
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=times(:,1,:)-1;
            days(pyr)=days(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            yr=times(:,1,:);
            days(nyr)=days(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        seconds=86400*mod(times(:,2,:),1);
    case 3
        % year, month, day of month
        ndays=[0 31 59 90 120 151 181 212 243 273 304 334].';
        times(:,1,:)=times(:,1,:)+floor((times(:,2,:)-1)/12);
        times(:,2,:)=mod(times(:,2,:)-1,12)+1;
        yr=times(:,1,:);
        days=yr*365+ndays(times(:,2,:))+floor(times(:,3,:))...
            +((~mod(yr,4)&(mod(yr,100)|~mod(yr,400)))&times(:,2,:)>2);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=yr-1;
            days(pyr)=days(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            days(nyr)=days(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        seconds=86400*mod(times(:,3,:),1);
    case 5
        % year, day of year, hour, minute, second
        days=times(:,1,:)*365+floor(times(:,2,:));
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=times(:,1,:)-1;
            days(pyr)=days(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            yr=times(:,1,:);
            days(nyr)=days(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        seconds=86400*mod(times(:,2,:),1)+times(:,3,:)*3600 ...
            +times(:,4,:)*60+times(:,5,:);
        days=days+floor(seconds/86400);
        seconds=mod(seconds,86400);
    case 6
        % year, month, day of month, hour, minute, second
        ndays=[0 31 59 90 120 151 181 212 243 273 304 334].';
        times(:,1,:)=times(:,1,:)+floor((times(:,2,:)-1)/12);
        times(:,2,:)=mod(times(:,2,:)-1,12)+1;
        yr=times(:,1,:);
        days=yr*365+ndays(times(:,2,:))+floor(times(:,3,:))...
            +((~mod(yr,4)&(mod(yr,100)|~mod(yr,400)))&times(:,2,:)>2);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=yr-1;
            days(pyr)=days(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            days(nyr)=days(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        seconds=86400*mod(times(:,3,:),1)+times(:,4,:)*3600 ...
            +times(:,5,:)*60+times(:,6,:);
        days=days+floor(seconds/86400);
        seconds=mod(seconds,86400);
end

% fix output to match
modserial=nan([sz(1) 2 sz(3:end)]);
modserial(:,1,:)=days; modserial(:,2,:)=seconds;

end
