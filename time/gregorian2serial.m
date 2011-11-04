function [serial]=gregorian2serial(times)
%GREGORIAN2SERIAL    Convert Gregorian dates to serial dates
%
%    Usage:    serialdates=gregorian2serial(gregoriandates)
%
%    Description:
%     GREGORIAN2SERIAL(DATES) returns the equivalent serial dates for the
%     Gregorian dates stored in DATES.  Serial date 1 corresponds to
%     January 1, 0000.  DATES must be a Nx2, Nx3, Nx5 or Nx6 numeric array
%     of either [year dayofyear], [year month dayofmonth],
%     [year dayofyear hour minute seconds] or
%     [year month dayofmonth hour minute seconds].
%
%    Notes:
%     - Does not account for UTC leap seconds
%     - Basically like Matlab's DATENUM except that it also accepts day of
%       year style dates, rolls-over negative calendar months, and does not
%       handle string input.  Because this is not compiled and DATENUM is,
%       DATENUM is 2 to 3 times faster.
%
%    Examples:
%     % The output of Matlab's datenum and gregorian2serial should match:
%     gregorian2serial([2008 366 0 0 0.01])
%     datenum([2008 12 31 0 0 0.01])
%
%    See also: SERIAL2GREGORIAN, GREGORIAN2MODSERIAL, MODSERIAL2GREGORIAN,
%              ISLEAPYEAR, DOY2CAL, CAL2DOY, FIXDATES, FIXTIMES, TIMEDIFF

%     Version History:
%        Nov.  2, 2008 - initial version
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
    error('seizmo:gregorian2serial:badInput',...
        'DATES must be a numeric Nx2, Nx3, Nx5 or Nx6 array!');
end

% work by date type
switch sz(2)
    case 2
        % year, day of year
        serial=times(:,1,:)*365+times(:,2,:);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=times(:,1,:)-1;
            serial(pyr)=serial(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            yr=times(:,1,:);
            serial(nyr)=serial(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
    case 3
        % year, month, day of month
        ndays=[0 31 59 90 120 151 181 212 243 273 304 334].';
        times(:,1,:)=times(:,1,:)+floor((times(:,2,:)-1)/12);
        times(:,2,:)=mod(times(:,2,:)-1,12)+1;
        yr=times(:,1,:);
        serial=yr*365+ndays(times(:,2,:))+times(:,3,:)...
            +((~mod(yr,4)&(mod(yr,100)|~mod(yr,400)))&times(:,2,:)>2);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=yr-1;
            serial(pyr)=serial(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            serial(nyr)=serial(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
    case 5
        % year, day of year, hour, minute, second
        serial=times(:,1,:)*365+times(:,2,:);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=times(:,1,:)-1;
            serial(pyr)=serial(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            yr=times(:,1,:);
            serial(nyr)=serial(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        serial=serial+times(:,3,:)/24+times(:,4,:)/1440+times(:,5,:)/86400;
    case 6
        % year, month, day of month, hour, minute, second
        ndays=[0 31 59 90 120 151 181 212 243 273 304 334].';
        times(:,1,:)=times(:,1,:)+floor((times(:,2,:)-1)/12);
        times(:,2,:)=mod(times(:,2,:)-1,12)+1;
        yr=times(:,1,:);
        serial=yr*365+ndays(times(:,2,:))+times(:,3,:)...
            +((~mod(yr,4)&(mod(yr,100)|~mod(yr,400)))&times(:,2,:)>2);
        pyr=times(:,1,:)>0;
        if(any(pyr))
            yr1=yr-1;
            serial(pyr)=serial(pyr)...
                +fix(yr1(pyr)/4)-fix(yr1(pyr)/100)+fix(yr1(pyr)/400)+1;
        end
        nyr=times(:,1,:)<0;
        if(any(nyr))
            serial(nyr)=serial(nyr)...
                +fix(yr(nyr)/4)-fix(yr(nyr)/100)+fix(yr(nyr)/400);
        end
        serial=serial+times(:,4,:)/24+times(:,5,:)/1440+times(:,6,:)/86400;
end

% fix output to match
serial=reshape(serial,[sz(1) 1 sz(3:end)]);

end
