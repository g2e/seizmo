function [leaps]=leapsinday(dates,option)
%LEAPSINDAY    Returns the number of leap seconds in the given dates
%
%    Usage:    seconds=leapsinday(dates)
%              seconds=leapsinday(dates,'serial')
%
%    Description: LEAPSINDAY(DATES) returns an array of additional seconds
%     to be added to each of the UTC dates in DATES.  If a positive leap
%     second occurs during one of the dates, that date will have a
%     corresponding 1 returned.  If there are no leap seconds during a
%     particular day, a 0 will be returned.  DATES must be a Nx2 array of
%     [years dayofyears] or a Nx3 array of [years months dayofmonths].
%
%     LEAPSINDAY(DATES,'SERIAL') takes a serialdate array as input.  Note
%     that this option really is not a good idea as using Matlab's DATENUM
%     will not properly handle UTC times (no leap second support) when
%     converting to a serial date number.  If your time is within a
%     positive leap second, DATENUM (and DATEVEC, DATESTR) will convert the
%     date to the following date thus making the date incorrect.
%
%    Notes:
%     - Only valid for UTC dates from 1972 on when the UTC second was
%       synced with the International Atomic Time (TAI) second and leap
%       seconds were introduced to keep UTC near UT1.
%
%    Examples:
%     Find out how many seconds there will be today:
%      86400+leapsinday(now,'serial')
%
%    See also: totalleaps, getleapseconds

%     Version History:
%        Nov.  1, 2008 - initial version
%        Nov. 10, 2008 - uses GETLEAPSECONDS to speed calls up
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - fix leap second bug, fix hardcode bug,
%                        minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:05 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end;

% check option
if(nargin==1 || isempty(option))
    option='gregorian';
elseif(~ischar(option) || ~any(strcmpi(option,{'serial' 'gregorian'})))
    error('seizmo:leapsinday:optionBad',...
        'OPTION must be ''serial'' or ''gregorian''!');
end

% require numeric
if(~isnumeric(dates))
    error('seizmo:leapsinday:badInput',...
        'DATES must be a numeric array!');
end

% work by input type
sz=size(dates);
switch lower(option)
    case 'gregorian'
        % require Nx2 or Nx3
        ndates=prod(sz)/sz(2);
        if(~any(sz(2)==[2 3]))
            error('seizmo:leapsinday:badInput',...
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
[leapdates,offset]=getleapseconds();
nleaps=numel(offset);

% correct for getleapseconds giving the end date of each leap second
leapdates=leapdates-1;

% expand (to vectorize)
dates=dates(:,ones(1,nleaps));
leapdates=leapdates(:,ones(1,ndates));

% compare dates
leaps=zeros(sz);
idx=(dates==leapdates.')*(1:nleaps).';
leaps(idx>0)=offset(idx(idx>0));

end
