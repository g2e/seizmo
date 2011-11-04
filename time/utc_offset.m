function [offsets]=utc_offset(dates,option)
%UTC_OFFSET    Returns the TAI-UTC offset for the given dates
%
%    Usage:    offsets=utc_offset(dates)
%              offsets=utc_offset(dates,'serial')
%
%    Description:
%     OFFSETS=UTC_OFFSET(DATES) returns the total offset in seconds of
%     International Atomic Time (TAI) from UTC time for any UTC date in the
%     array OFFSETS. DATES must be an array of UTC dates of the form
%     [years daysofyear] or [years months daysofmonth].
%
%     OFFSETS=UTC_OFFSET(DATES,'SERIAL') takes a serialdate array as input.
%     Note that this option really is not a good idea as using Matlab's
%     DATENUM will not properly handle UTC times (no leap second support)
%     when converting to a serial date number.  If your time is within a
%     positive leap second, DATENUM (and DATEVEC, DATESTR) will convert the
%     date to the following date thus making the date incorrect.
%
%    Notes:
%     - Offsets are not valid for pre-1972 when UTC was not synced with the
%       TAI second and drift between the two time standards occured. An
%       error message now occurs for pre-1972 dates.
%
%    Examples:
%     % Get the current TAI-UTC offset:
%     utc_offset(now,'serial')
%
%     % The GPS-UTC offset is a slight modification (always 19s less):
%     utc_offset(now,'serial')-19
%
%    See also: UTC_LOD, LEAPSECONDS, LEAPSECONDS_UPDATE, FIXTIMES,
%              TIMEDIFF, UTC2TAI, TAI2UTC, ISLEAPYEAR, FIXDATES

%     Version History:
%        Nov.  1, 2008 - initial version
%        Nov. 10, 2008 - uses GETLEAPSECONDS to speed calls up
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - minor doc update
%        Sep.  5, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - rename from TOTALLEAPS to UTC_OFFSET, added error
%                        message if pre-1972 dates are found, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check option
if(nargin==1 || isempty(option))
    option='gregorian';
elseif(~ischar(option) || ~isvector(option) ...
        || ~any(strcmpi(option,{'serial' 'gregorian'})))
    error('seizmo:utc_offset:optionBad',...
        'OPTION must be ''serial'' or ''gregorian''!');
end

% require numeric
if(~isnumeric(dates))
    error('seizmo:utc_offset:badInput',...
        'DATES must be a numeric array!');
end

% work by input type
sz=size(dates);
switch lower(option)
    case 'gregorian'
        % require Nx2 or Nx3
        ndates=prod(sz)/sz(2);
        if(~any(sz(2)==[2 3]))
            error('seizmo:utc_offset:badInput',...
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

% check for pre-1972 dates
if(any(dates<datenum([1972 1 1])))
    error('seizmo:utc_offset:unsupportedDate',...
        'Pre-1972 date found: UTC dates prior to 1972 are unsupported!');
end

% get leapsecond info
[leapdates,leaps,offsets]=leapseconds;
nleaps=numel(leaps);
offsets=[nan; offsets]; % dates prior have undefined offset currently

% expand (to vectorize)
dates=dates(:,ones(1,nleaps));
leapdates=leapdates(:,ones(1,ndates));

% compare dates
offsets=reshape(offsets(1+sum(dates>=leapdates.',2)),sz);

end
