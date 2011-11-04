function [seconds]=utc_lod(dates,option)
%UTC_LOD    Returns the number of seconds in the given days
%
%    Usage:    seconds=utc_lod(dates)
%              seconds=utc_lod(dates,'serial')
%
%    Description:
%     SECONDS=UTC_LOD(DATES) returns an array of the number of seconds for
%     each of the UTC dates in DATES.  If a positive leap second occurs
%     during one of the dates, that date will have a corresponding 86401
%     returned.  If there are no leap seconds during a particular day,
%     86400 will be returned.  DATES must be a Nx2 array of
%     [years dayofyears] or a Nx3 array of [years months dayofmonths].
%
%     SECONDS=UTC_LOD(DATES,'SERIAL') takes a serialdate array as input.
%     Note that this option really is not a good idea as using Matlab's
%     DATENUM will not properly handle UTC times (no leap second support)
%     when converting to a serial date number.  If your time is within a
%     positive leap second, DATENUM (and DATEVEC, DATESTR) will convert the
%     date to the following date thus making the date incorrect.
%
%    Notes:
%     - Length of day results are not valid for pre-1972 when UTC was not synced with the
%       TAI second and drift between the two time standards occured. An
%       error message now occurs for pre-1972 dates.
%
%    Examples:
%     % Find out how many seconds there will be today (as
%     % long as you are not currently in the leap second):
%     utc_lod(now,'serial')
%
%    See also: UTC_OFFSET, LEAPSECONDS, LEAPSECONDS_UPDATE, FIXTIMES,
%              TIMEDIFF, UTC2TAI, TAI2UTC, ISLEAPYEAR, FIXDATES

%     Version History:
%        Nov.  1, 2008 - initial version
%        Nov. 10, 2008 - uses GETLEAPSECONDS to speed calls up
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - fix leap second bug, fix hardcode bug,
%                        minor doc update
%        Sep.  5, 2009 - minor doc update
%        Oct. 16, 2009 - major speed increase
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - rename from LEAPSINDAY to UTC_LOD, added error
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
elseif(~ischar(option) || ~any(strcmpi(option,{'serial' 'gregorian'})))
    error('seizmo:utc_lod:optionBad',...
        'OPTION must be ''serial'' or ''gregorian''!');
end

% require numeric
if(~isnumeric(dates))
    error('seizmo:utc_lod:badInput',...
        'DATES must be a numeric array!');
end

% work by input type
sz=size(dates);
switch lower(option)
    case 'gregorian'
        % require Nx2 or Nx3
        ndates=prod(sz)/sz(2);
        if(~any(sz(2)==[2 3]))
            error('seizmo:utc_lod:badInput',...
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
    error('seizmo:utc_lod:unsupportedDate',...
        'Pre-1972 date found: UTC dates prior to 1972 are unsupported!');
end

% get leapsecond info
[leapdates,leaps]=leapseconds;
nleaps=numel(leaps);

% correct for leapseconds giving the end date of each leap second
leapdates=leapdates-1;

% expand (to vectorize)
dates=dates(:,ones(1,nleaps));
leapdates=leapdates(:,ones(1,ndates));

% compare dates
seconds=86400*ones(sz);
[row,col]=find(dates==leapdates.');
if(~isempty(row)); seconds(row)=seconds(row)+leaps(col); end

end
