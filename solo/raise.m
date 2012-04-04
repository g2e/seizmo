function [data]=raise(data,power)
%RAISE    Scales the dependent component of SEIZMO records by a given power
%
%    Usage:    data=raise(data,power)
%
%    Description:
%     DATA=RAISE(DATA,POWER) raises SEIZMO records in DATA by POWER.  POWER
%     must be a real-valued scalar such as 3 or 1/4.  The default value for
%     POWER is 1 which leaves the records unchanged.  Note that the sign of
%     the data is preserved and the power just scales the absolute value of
%     the data.
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
%     % Cube the data in records:
%     data=raise(data,3);
%
%     % Get the 4th root of records:
%     data=raise(data,1/4);
%
%    See also: SOLOFUN

%     Version History:
%        Mar. 16, 2010 - initial version
%        Jan.  6, 2011 - drop versioninfo caching, fix nargchk, update for
%                        seizmofun/solofun name change
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 17:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt raising record's to power
try
    % default power
    if(isempty(power)); power=1; end
    
    % check power
    if(~isscalar(power) || ~isreal(power))
        error('seizmo:raise:badInput',...
            'RAISE must be a real-valued scalar!');
    end
    
    % send to solofun
    data=solofun(data,@(x)sign(x).*abs(x).^power);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
