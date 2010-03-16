function [data]=raise(data,power)
%RAISE    Scales the dependent component of SEIZMO records by a given power
%
%    Usage:    data=raise(data,power)
%
%    Description: DATA=RAISE(DATA,POWER) raises SEIZMO records in DATA by
%     POWER.  POWER must be a real-valued scalar such as 3 or 1/4.  The
%     default value for POWER is 1 which leaves the records unchanged.
%     Note that the sign of the data is preserved and the power just scales
%     the absolute value of the data.
%
%    Notes:
%     - POWER > 1 emphasizes strong peaks
%     - 0 < POWER < 1 emphasizes smaller peaks
%     - POWER==0 is a 1-bit filter
%     - POWER < 0 is a bad idea!
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     Cube the data in records:
%      data=raise(data,3);
%
%     Get the 4th root of records:
%      data=raise(data,1/4);
%
%    See also: SEIZMOFUN

%     Version History:
%        Mar. 16, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 16, 2010 at 15:45 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt raising record's to power
try
    % default power
    if(isempty(power)); power=1; end
    
    % check power
    if(~isscalar(power) || ~isreal(power))
        error('seizmo:raise:badInput',...
            'RAISE must be a real-valued scalar!');
    end
    
    % send to seizmofun
    data=seizmofun(data,@(x)sign(x).*abs(x).^power);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
