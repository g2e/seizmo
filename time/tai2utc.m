function [tai]=tai2utc(tai)
%TAI2UTC    Convert TAI time to UTC time
%
%    Usage:    utc=tai2utc(tai)
%
%    Description:
%     UTC=TAI2UTC(TAI) returns the equivalent Universal Coordinated Times
%     (UTC) of the International Atomic Times (TAI) in TAI.  Times should
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
%     % TAI to UTC can handle leap seconds:
%     tai2utc([2009 1 1 0 0 22; 2009 1 1 0 0 23; 2009 1 1 0 0 24])
%
%    See also: UTC2TAI, FIXTIMES, TIMEDIFF, LEAPSECONDS, UTC_OFFSET,
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
sz=size(tai);
if(~isnumeric(tai) || ~any(sz(2)==[5 6]))
    error('seizmo:tai2utc:badInput','TAI array must be Nx5 or Nx6!');
end

% fix times
tai=fixtimes(tai);

% UTC time
tai(:,end,:)=tai(:,end,:)-utc_offset(tai(:,1:end-3,:));

% fix times
tai=fixtimes(tai,'utc');

end
