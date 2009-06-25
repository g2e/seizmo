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
%    See also: getheader, changeheader, listheader

%     Version History:
%        June 25, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 25, 2009 at 10:30 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

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
    % check headers
    data=checkheader(data);
    
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
set_seizmocheck_state(oldseizmocheckstate);

end
