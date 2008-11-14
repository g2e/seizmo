function [times]=fixtimes(times,option)
%FIXTIMES    Cleans up times to be in the proper ranges
%
%    Description: FIXTIMES(TIMES) returns equivalent times to the times in
%     TIMES such that they follow typical Gregorian calendar conventions
%     (months are 1 to 12, dayofmonths are 1 to the last day of the month,
%     and dayofyears are from 1 to 365 or 366 if it is a leap year) and
%     the 24 hour style of time (hours 0-23, minutes 0-59 and seconds
%     0-60).  TIMES must be a Nx5 array of [yr dayofyr hr min sec] or a Nx6
%     array of [yr mon dayofmon hr min sec].  Only the seconds portion
%     TIMES is allowed to be non-integer (ie you cannot have 1.5 minutes
%     etc).
%
%     FIXTIMES(TIMES,'UTC') allows fixing UTC times which have leap seconds
%     occasionally inserted on certain dates (see LEAPSECONDS).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    times=fixtimes(times)
%              times=fixtimes(times,'utc')
%
%    Examples:
%     Find the UTC time 500 seconds after some other time in UTC:
%      UTC=[2008 12 31 23 55 3.789];
%      UTC(6)=UTC(6)+500;
%      fixtimes(UTC,'UTC')
%
%    See also: fixdates, isleapyear, leapseconds, totalleaps, leapsinday

%     Version History:
%        Nov.  2, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  2, 2008 at 08:20 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check times
sz=size(times);
if(~isnumeric(times) || ~any(sz(2)==[5 6]))
    error('SAClab:fixtimes:badInput',...
        'TIMES must be a numeric array!');
end

% check option
if(nargin==1 || isempty(option))
    option='tai';
elseif(~ischar(option) || ~any(strcmpi(option,{'utc' 'tai'})))
    error('SAClab:fixtimes:optionBad',...
        'OPTION must be ''utc'' or ''tai''!');
end

% permute and expand
times=permute(times,[2 1 3:numel(sz)]);
times=times(:,:).';

% proceed by option
switch lower(option)
    case 'tai'
        % fix TAI time
        times(:,end-1)=times(:,end-1)+floor(times(:,end)/60);
        times(:,end)=mod(times(:,end),60);
        times(:,end-2)=times(:,end-2)+floor(times(:,end-1)/60);
        times(:,end-1)=mod(times(:,end-1),60);
        times(:,end-3)=times(:,end-3)+floor(times(:,end-2)/24);
        times(:,end-2)=mod(times(:,end-2),24);
        times(:,1:end-3)=fixdates(times(:,1:end-3));
    case 'utc'
        % fix UTC time (to the minute)
        times(:,end-2)=times(:,end-2)+floor(times(:,end-1)/60);
        times(:,end-1)=mod(times(:,end-1),60);
        times(:,end-3)=times(:,end-3)+floor(times(:,end-2)/24);
        times(:,end-2)=mod(times(:,end-2),24);
        times(:,1:end-3)=fixdates(times(:,1:end-3));
        
        % get what second we are at in the day
        secnow=times(:,end)+times(:,end-1)*60+times(:,end-2)*3600;
        sectoday=86400+leapsinday(times(:,1:end-3));
        
        % rollover seconds if needed
        low=find(secnow<0);
        while(~isempty(low))
            times(low,end-3)=times(low,end-3)-1;
            sectoday(low)=86400+leapsinday(times(low,1:end-3));
            secnow(low)=secnow(low)+sectoday(low);
            low=low(secnow<0);
        end
        high=find(secnow>=sectoday);
        while(~isempty(high))
            times(high,end-3)=times(high,end-3)+1;
            secnow(high)=secnow(high)-sectoday(high);
            sectoday(high)=86400+leapsinday(times(high,1:end-3));
            high=high(secnow(high)>=sectoday(high));
        end
        
        % clean up times one last time
        times(:,1:end-3)=fixdates(times(:,1:end-3));
        times(:,end-2)=floor(secnow/3600);
        times(times(:,end-2)>23,end-2)=23;
        secnow=secnow-times(:,end-2)*3600;
        times(:,end-1)=floor(secnow/60);
        times(times(:,end-1)>59,end-1)=59;
        times(:,end)=secnow-times(:,end-1)*60;
end

% reshape back
times=permute(reshape(times.',sz([2 1 3:end])),[2 1 3:numel(sz)]);

end
