function [leaps]=leapsinday(dates,option)
%LEAPSINDAY    Returns the number of leap seconds in the given dates
%
%    Description: LEAPSINDAY(DATES) returns the an array of additional
%     seconds added to each of the UTC dates in DATES.  If a positive leap
%     second occurs during one of the dates, that date will have a
%     corresponding 1 returned.  If there are no leap seconds during a
%     particular day, a 0 will be returned.  DATES must be a Nx2 array of
%     [years dayofyears] or a Nx3 array of [years months dayofmonths].
%
%     LEAPSINDAY(DATES,'SERIAL') takes a serialdate array as input.  Note
%     that this option really is not a good idea as using Matlab's DATENUM
%     will not properly handle UTC times (no leap second support) when
%     converting to a serial date number.
%
%    Notes:
%     - Only valid for UTC dates from 1972 on when the UTC second was
%       synced with the International Atomic Time (TAI) second and leap
%       seconds were introduced to keep UTC near UT1.
%
%    Tested on: Matlab r2007b
%
%    Usage:    seconds=leapsinday(dates)
%              seconds=leapsinday(dates,'serial')
%
%    Examples:
%     Find out how many seconds there will be in today:
%      86400+leapsinday(now,'serial')
%
%    See also: totalleaps, leapseconds

%     Version History:
%        Nov.  1, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2008 at 17:30 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check option
if(nargin==1 || isempty(option))
    option='gregorian';
elseif(~ischar(option) || ~any(strcmpi(option,{'serial' 'gregorian'})))
    error('SAClab:leapsinday:optionBad',...
        'OPTION must be ''serial'' or ''gregorian''!');
end

% require numeric
if(~isnumeric(dates))
    error('SAClab:leapsinday:badInput',...
        'DATES must be a numeric array!');
end

% work by input type
sz=size(dates);
switch lower(option)
    case 'gregorian'
        % require Nx2 or Nx3
        ndates=prod(sz)/sz(2);
        if(~any(sz(2)==[2 3]))
            error('SAClab:leapsinday:badInput',...
                'DATES must be a Nx2 or Nx3 array!');
        end
        
        % clean up input
        dates=fixdates(dates);
        if(sz(2)==2); dates=doy2cal(dates); end
        
        % permute and pass to datenum
        dates=permute(dates,[2 1 3:numel(sz)]);
        dates=datenum(dates(:,:).');
        sz(2)=1;
    case 'serial'
        % integer only
        ndates=prod(sz);
        dates=floor(dates(:));
end
        
% get leapsecond info
leapstr=leapseconds;
nleaps=size(leapstr,1);
leapdates=datenum(leapstr)-1;
offset=2*strcmp('+',cellstr(leapstr(:,end)))-1;

% expand (to vectorize)
dates=dates(:,ones(1,nleaps));
leapdates=leapdates(:,ones(1,ndates));

% compare dates
leaps=zeros(sz);
idx=(dates==leapdates.')*(1:24).';
leaps(idx>0)=offset(idx(idx>0));

end
