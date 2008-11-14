function [diff]=timediff(times1,times2,option)
%TIMEDIFF    Return number of seconds between times
%
%    Description: TIMEDIFF(TIMES1,TIMES2) returns the difference in time
%     between times in TIMES1 and times in TIMES2.  TIMES1 and TIMES2 must
%     be a Nx2, Nx3, Nx5, or Nx6 array of [yr dayofyr], [yr mon dayofmon], 
%     [yr dayofyr hr min sec] or [yr mon dayofmon hr min sec].  Only the
%     seconds portion of TIMES1 and TIMES2 is allowed to be non-integer
%     (ie you cannot have 1.5 minutes etc).  TIMES1 and TIMES2 must have
%     equal size of be a single time.  Time difference is returned in
%     number of seconds.
%
%     TIMEDIFF(TIMES1,TIMES2,'UTC') allows finding the difference between
%     UTC times (which may have leap seconds occasionally inserted on
%     certain dates -- see LEAPSECONDS).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    diff=timediff(times1,times2)
%              diff=timediff(times1,times2,'utc')
%
%    Examples:
%     Find the number of seconds in 2005:
%      timediff([2005 1],[2006 1],'utc')
%
%    See also: gregorian2modserial, utc2tai

%     Version History:
%        Nov. 12, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 12, 2008 at 04:40 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin))

% check times
sz1=size(times1);
sz2=size(times2);
if(~isnumeric(times1) || ~isnumeric(times2)...
        || isempty(times1) || isempty(times2)...
        || ~any(sz1(2)==[2 3 5 6]) || ~any(sz2(2)==[2 3 5 6])...
        || (~isequal(sz1,sz2) && ...
        prod(sz1([1 3:end]))~=1 && prod(sz2([1 3:end]))~=1))
    error('SAClab:timediff:badInput',...
        'TIMES1 and TIMES2 must be a non-empty numeric date arrays!');
end

% check option
if(nargin==2 || isempty(option))
    option='tai';
elseif(~ischar(option) || ~any(strcmpi(option,{'utc' 'tai'})))
    error('SAClab:timediff:optionBad',...
        'OPTION must be ''utc'' or ''tai''!');
end

% proceed by option
switch lower(option)
    case 'tai'
        % TAI => MODIFIED SERIAL
        modserial=gregorian2modserial(times2)-gregorian2modserial(times1);
        
        % take difference
        diff=submat(modserial,2,1)*86400+submat(modserial,2,2);
    case 'utc'
        % add zeros (date => time)
        if(any(sz1(2)==[2 3])); times1(:,end+3,:)=0; end
        if(any(sz2(2)==[2 3])); times2(:,end+3,:)=0; end
        
        % UTC => TAI => MODIFIED SERIAL
        modserial=gregorian2modserial(utc2tai(times2))...
            -gregorian2modserial(utc2tai(times1));
        
        % take difference
        diff=submat(modserial,2,1)*86400+submat(modserial,2,2);
end

end
