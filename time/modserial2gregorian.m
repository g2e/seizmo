function [times]=modserial2gregorian(modserial,option)
%MODSERIAL2GREGORIAN    Convert modified serial dates to Gregorian dates
%
%    Usage:    gregoriandates=modserial2gregorian(modserialdates)
%              gregoriandates=modserial2gregorian(modserialdates,option)
%
%    Description:
%     MODSERIAL2GREGORIAN(DATES) returns the equivalent Gregorian dates in
%     [year month dayofmonth hour minute second] format for the modified
%     serial dates stored in MODSERIAL.  Modified serial dates are the
%     number of days and seconds since January 0, year 0 at 00:00:00 and
%     are in the format [days seconds]. This is useful for maintaining
%     accuracy on subsecond operations.
%
%     MODSERIAL2GREGORIAN(DATES,OPTION) specifies the output format:
%      OPTION        OUTPUT
%       'caldate' =>  [year month dayofmonth]
%       'caltime' =>  [year month dayofmonth hour minute seconds]
%       'doydate' =>  [year dayofyear]
%       'doytime' =>  [year dayofyear hour minute seconds]
%
%    Notes:
%     - Does not account for UTC leap seconds
%     - Basically like Matlab's DATEVEC except that it takes in modified
%       serial dates and has an option to output day-of-year style dates.
%       It also does not handle string input and isn't compiled so it is
%       a bit slower.
%
%    Examples:
%     % 500 seconds from now:
%     modserial2gregorian([now 500])
%
%    See also: GREGORIAN2MODSERIAL, SERIAL2GREGORIAN, GREGORIAN2SERIAL,
%              FIXTIMES, FIXDATES, TIMEDIFF, ISLEAPYEAR, CAL2DOY, DOY2CAL

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Sep. 23, 2009 - fixed serial conversion (year 0 bug)
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - doc update
%        Feb. 13, 2012 - vector input for cal* output bugfix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check days and seconds
sz=size(modserial);
if(~isnumeric(modserial) || sz(2)~=2)
    error('seizmo:modserial2gregorian:badInput',...
        'MODSERIAL must be a numeric array of size Nx2!');
end

% check option
if(nargin==1 || isempty(option))
    option='caltime';
elseif(~ischar(option)...
        || ~any(strcmpi(option,{'caldate' 'caltime' 'doydate' 'doytime'})))
    error('seizmo:modserial2gregorian:optionBad',...
        ['OPTION must be ''caldate'', ''caltime'', '...
        '''doydate'' or ''doytime''!']);
end

% clean up
modserial(:,2,:)=modserial(:,2,:)+mod(modserial(:,1,:),1)*86400;
modserial(:,1,:)=floor(modserial(:,1,:))+floor(modserial(:,2,:)/86400);
modserial(:,2,:)=mod(modserial(:,2,:),86400);

% get year
serial=modserial(:,1,:);
yr=floor(serial/365.2425);
newserial=serial-(365*yr+fix((yr+3)/4)-fix((yr+99)/100)+fix((yr+399)/400));
low=newserial<1;
if(any(low(:)))
    yr(low)=yr(low)-1;
    serial(low)=serial(low)-(365*yr(low)+fix((yr(low)+3)/4) ...
        -fix((yr(low)+99)/100)+fix((yr(low)+399)/400));
end
serial(~low)=newserial(~low);
clear newserial

% convert
switch lower(option)
    case 'caldate'
        times=nan([sz(1) 3 sz(3:end)]);
        ndays=[1 32 60 91 121 152 182 213 244 274 305 335];
        times(:,1,:)=yr;
        ndays2=ndays(ones(sz(1),1),:,ones(prod(sz(3:end)),1));
        leap=isleapyear(times(:,1,:));
        ndays2(:,3:end,:)=ndays2(:,3:end,:)+leap(:,ones(1,10),:);
        times(:,2,:)=sum(serial(:,ones(12,1),:)>=ndays2,2);
        times(:,3,:)=serial+modserial(:,2,:)/86400 ...
            -reshape(ndays(times(:,2,:)),size(times(:,2,:)))...
            -(leap & times(:,2,:)>2)+1;
    case 'caltime'
        times=nan([sz(1) 6 sz(3:end)]);
        ndays=[1 32 60 91 121 152 182 213 244 274 305 335];
        times(:,1,:)=yr;
        ndays2=ndays(ones(sz(1),1),:,ones(prod(sz(3:end)),1));
        leap=isleapyear(times(:,1,:));
        ndays2(:,3:end,:)=ndays2(:,3:end,:)+leap(:,ones(1,10),:);
        times(:,2,:)=sum(serial(:,ones(12,1),:)>=ndays2,2);
        times(:,3,:)=serial-reshape(ndays(times(:,2,:)),...
            size(times(:,2,:)))-(leap & times(:,2,:)>2)+1;
        times(:,4,:)=fix(modserial(:,2,:)/3600);
        modserial(:,2,:)=modserial(:,2,:)-times(:,4,:)*3600;
        times(:,5,:)=fix(modserial(:,2,:)/60);
        modserial(:,2,:)=modserial(:,2,:)-times(:,5,:)*60;
        times(:,6,:)=modserial(:,2,:);
    case 'doydate'
        times=nan([sz(1) 2 sz(3:end)]);
        times(:,1,:)=yr;
        times(:,2,:)=serial+modserial(:,2,:)/86400;
    case 'doytime'
        times=nan([sz(1) 5 sz(3:end)]);
        times(:,1,:)=yr;
        times(:,2,:)=serial;
        times(:,3,:)=fix(modserial(:,2,:)/3600);
        modserial(:,2,:)=modserial(:,2,:)-times(:,3,:)*3600;
        times(:,4,:)=fix(modserial(:,2,:)/60);
        modserial(:,2,:)=modserial(:,2,:)-times(:,4,:)*60;
        times(:,5,:)=modserial(:,2,:);
end

end
