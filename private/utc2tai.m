function [utc]=utc2tai(utc)
%UTC2TAI    Convert UTC time to TAI time
%
%    Description: UTC2TAI(UTC) returns the equivalent International Atomic
%     Time (TAI) of the UTC times in UTC.  Times should be either a Nx5 or
%     Nx6 of [yr dayofyr hr min sec] or [yr mon dayofmon hr min sec].
%
%    Notes:
%     - Only valid for UTC dates from 1972 on when the UTC second was
%       synced with the International Atomic Time (TAI) second and leap
%       seconds were introduced to keep UTC near UT1.
%
%    Tested on: Matlab r2007b
%
%    Usage:    tai=utc2tai(utc)
%
%    Examples:
%     It is a lot easier to do time differences in TAI,
%     which does not have leap seconds:
%      utc2tai([2008 12 31 23 59 60])-utc2tai([2009 1 1 0 0 0])
%
%    See also: tai2utc, fixtimes, timediff, leapseconds, totalleaps

%     Version History:
%        Nov.  2, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  2, 2008 at 15:10 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require Nx5 or Nx6 numeric
sz=size(utc);
if(~isnumeric(utc) || ~any(sz(2)==[5 6]))
    error('SAClab:utc2tai:badInput','UTC array must be Nx5 or Nx6!');
end

% fix times
utc=fixtimes(utc,'utc');

% TAI time
utc(:,end,:)=utc(:,end,:)+totalleaps(utc(:,1:end-3,:));

% fix times
utc=fixtimes(utc);

end
