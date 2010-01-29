function [data]=changetype(data,type)
%CHANGETYPE    Changes file type of SEIZMO records
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Header changes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Dec. 30, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 30, 2009 at 04:50 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt type change
try
    % check type
    
    % number of records
    nrecs=numel(data);
    
    % loop over records
    for i=1:nrecs
        % act by type
        switch type
            case 'sac'
            case 'seizmo'
        end
    end
catch
    
end

% sac to seizmo
% - allow ncmp
% - double precision

% seizmo to sac
% - drop ncmp
% - single precision (skip)

end