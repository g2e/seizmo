function [data]=copyheader(data,idx,varargin)
%COPYHEADER    Copy one record's header to all records
%
%    Usage:    data=copyheader(data)
%              data=copyheader(data,idx)
%              data=copyheader(data,name)
%              data=copyheader(data,idx|name,field1,...,fieldN)
%
%    Description: COPYHEADER(DATA) copies the entire header of record 1 in
%     DATA to the rest of the records in DATA.
%
%     COPYHEADER(DATA,IDX) copies the entire header of the record indicated
%     by IDX to the rest of the records in DATA.
%
%     COPYHEADER(DATA,NAME) allows choosing the record by name rather than
%     by index.  NAME is compared to the name field in SEIZMO struct DATA
%     and the first record with a matching name is chosen.
%
%     COPYHEADER(DATA,IDX|NAME,FIELD1,...,FIELDN) copies fields FIELD1 to
%     FIELDN from the indicated record to the rest of the records in DATA.
%
%    Notes:
%
%    Header changes: POTENTIALLY ALL
%
%    Examples:
%     Copy picks from record 2 to all records:
%      data=copyheader(data,2,'t','kt')
%
%     Copy event info from WUSTL.sac to all other records:
%      data=copyheader(data,'WUSTL.sac','evla','evlo','evel','evdp');
%
%    See also: GETHEADER, CHANGEHEADER, LISTHEADER

%     Version History:
%        June 25, 2009 - initial version
%        Jan. 28, 2010 - seizmoverbose support, proper SEIZMO handling
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2010 at 21:20 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt copy of header
try
    % detail message
    if(seizmoverbose)
        disp('Copying Header of One Record to the Rest');
    end
    
    % number of records
    nrecs=numel(data);

    if(nargin==1)
        [data.head]=deal(getheader(data(1)));
    elseif(nargin==2)
        % check index
        if(isempty(idx))
            idx=1;
        elseif(ischar(idx) && size(idx,1)==1)
            idx=find(strcmp(idx,{data.name}),1);
            if(isempty(idx))
                error('seizmo:copyheader:badName','NAME not found!');
            end
        end
        if(~isscalar(idx) || ~isnumeric(idx) ...
                || idx~=fix(idx) || idx<1 || idx>nrecs)
            error('seizmo:copyheader:badIDX',...
                'IDX must be a valid scalar indice!');
        end

        % copy full header
        [data.head]=deal(getheader(data(idx)));
    else
        % check index
        if(isempty(idx))
            idx=1;
        elseif(ischar(idx) && size(idx,1)==1)
            idx=find(strcmp(idx,{data.name}),1);
            if(isempty(idx))
                error('seizmo:copyheader:badName','NAME not found!');
            end
        end
        if(~isscalar(idx) || ~isnumeric(idx) ...
                || idx~=fix(idx) || idx<1 || idx>nrecs)
            error('seizmo:copyheader:badIDX',...
                'IDX must be a valid scalar indice!');
        end

        % get values
        values=cell(1,numel(varargin));
        [values{:}]=getheader(data(idx),varargin{:});

        % apply values
        values=[varargin; values];
        data=changeheader(data,values{:});
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
