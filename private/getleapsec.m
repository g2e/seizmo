function [leap]=getleapsec(minute)
%GETLEAPSEC    Returns the seconds to add/subtract for specific UTC minutes
%
%    Description: GETLEAPSEC(MINUTES) takes an Nx4 numeric array of
%     [year dayofyear hour minute] and returns a value for each row
%     indicating whether there is a positive leap second (1), a negative
%     leap second (-1) or no leap second (0) during that UTC minute.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    leaps=getleapsec(minutes)
%
%    Examples:
%     How many seconds are in a specific UTC minute:
%      60+getleapsec([2008 366 23 59])
%
%    See also: leapseconds, utc2utc

%     Version History:
%        Oct. 30, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 29, 2008 at 20:10 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% transpose minute
minute=permute(minute,[2 1 3:ndims(minute)]);

% check input
[sz1,sz2]=size(minute);
if(~isnumeric(minute) || sz1~=4)
    error('SAClab:getleapsec:badInput','Input must be [year doy hr min]!');
end

% check for nonintegers
if(any(minute~=fix(minute)))
    error('SAClab:getleapsec:badInput','Input must be integers!');
end

% check range
badjday=(minute(2,:)<1 | minute(2,:)>366 ...
    | (minute(2,:)==366 & ~isleapyear(minute(1,:))));
badhour=(minute(3,:)<0 | minute(3,:)>23);
badmin=(minute(4,:)<0 | minute(4,:)>59);
if(any(badjday | badhour | badmin))
    error('SAClab:getleapsec:badRange','Time not valid!');
end

% calday
cdate=calday(minute(1:2,:).');

% get serial date
serialminute=datenum([cdate minute(3:4,:).' zeros(sz2,1)]);

% get leapsecond info
leapstr=leapseconds;
leapminute=datenum(leapstr)-60/86400;
offset=2*strcmp('+',cellstr(leapstr(:,end)))-1;

% find matches
[i,j]=ismember(serialminute,leapminute);

% get leap
sz=size(minute);
leap=zeros([sz(2) 1 sz(3:end)]);
leap(i)=offset(j(i));

end
