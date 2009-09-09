function [utc]=utc2tai(utc)
%UTC2TAI    Convert UTC time to TAI time
%
%    Usage:    tai=utc2tai(utc)
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
%    Examples:
%     It is a lot easier to do time differences in TAI,
%     which does not have leap seconds:
%      utc2tai([2008 12 31 23 59 60])-utc2tai([2009 1 1 0 0 0])
%
%    See also: tai2utc, fixtimes, timediff, leapseconds, totalleaps,
%              getleapseconds, leapsinday

%     Version History:
%        Nov.  2, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2009 at 19:30 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end;

% require Nx5 or Nx6 numeric
sz=size(utc);
if(~isnumeric(utc) || ~any(sz(2)==[5 6]))
    error('seizmo:utc2tai:badInput','UTC array must be Nx5 or Nx6!');
end

% fix times
utc=fixtimes(utc,'utc');

% TAI time
utc(:,end,:)=utc(:,end,:)+totalleaps(utc(:,1:end-3,:));

% fix times
utc=fixtimes(utc);

end
