function [utc]=utc2tai(utc)
%UTC2TAI    Convert UTC time to TAI time
%
%    Usage:    tai=utc2tai(utc)
%
%    Description:
%     TAI=UTC2TAI(UTC) returns the equivalent International Atomic Time
%     (TAI) of the Universal Coordinated Times (UTC) in UTC.  Times should
%     be either a Nx5 or Nx6 of [yr dayofyr hr min sec] or
%     [yr mon dayofmon hr min sec].
%
%    Notes:
%     - Only valid for UTC dates from 1972 on when the UTC second was
%       synced with the International Atomic Time (TAI) second and leap
%       seconds were introduced to keep UTC near UT1. An error is issued
%       for pre-1972 dates.
%
%    Examples:
%     % It is a lot easier to do time differences in TAI,
%     % which does not have leap seconds:
%     utc2tai([2008 12 31 23 59 60])-utc2tai([2009 1 1 0 0 0])
%
%    See also: TAI2UTC, FIXTIMES, TIMEDIFF, LEAPSECONDS, UTC_OFFSET,
%              LEAPSECONDS_UPDATE, UTC_LOD, ISLEAPYEAR, FIXDATES

%     Version History:
%        Nov.  2, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - doc update, fix for UTC_OFFSET
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require Nx5 or Nx6 numeric
sz=size(utc);
if(~isnumeric(utc) || ~any(sz(2)==[5 6]))
    error('seizmo:utc2tai:badInput','UTC array must be Nx5 or Nx6!');
end

% fix times
utc=fixtimes(utc,'utc');

% TAI time
utc(:,end,:)=utc(:,end,:)+utc_offset(utc(:,1:end-3,:));

% fix times
utc=fixtimes(utc);

end
