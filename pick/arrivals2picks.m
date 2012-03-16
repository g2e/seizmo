function [data]=arrivals2picks(data,phaselist,idx,varargin)
%ARRIVALS2PICKS    Push arrivals to header pick fields in SEIZMO records
%
%    Usage:    data=arrivals2picks(data)
%              data=arrivals2picks(data,phaselist)
%              data=arrivals2picks(data,phaselist,idx)
%              data=arrivals2picks(data,phaselist,idx,k)
%              data=arrivals2picks(data,phaselist,idx,k,'last')
%
%    Description:
%     DATA=ARRIVALS2PICKS(DATA) inserts all arrivals into the pick header
%     fields (KT/T/USER) until there are no more arrivals or the header
%     fields are exhausted (this is the usual case as there are only 10).
%     KT fields are given the phase name, T fields contain the travel time
%     corrected for the origin time, & USER fields give the ray parameter
%     in seconds per degree.
%
%     DATA=ARRIVALS2PICKS(DATA,PHASELIST) inserts the specified arrivals
%     into the pick header fields (KT/T/USER).  PHASELIST may be either a
%     cell array of phases or a comma-separated string.
%
%     DATA=ARRIVALS2PICKS(DATA,PHASELIST,IDX) indicates the indices
%     of the KT/T/USER header fields to insert arrivals.  In other words,
%     these are the fields than are possible overwrote.  The default is
%     0:9.  The indices are for all records and can not be individually
%     set.
%
%     DATA=ARRIVALS2PICKS(DATA,PHASELIST,IDX,K) specifies the maximum
%     number of phases to insert for each record.  Inserts the first K
%     arrivals found.
%
%     DATA=ARRIVALS2PICKS(DATA,PHASELIST,IDX,K,'LAST') inserts the last K
%     arrivals found matching the PHASELIST.
%
%    Notes:
%
%    Examples:
%     % Add 1st P arrivals to the last field index:
%     data=arrivals2picks(data,'P',9,1);
%
%    See also: MAKEARRIVALS, FINDARRIVALS, FINDPICKS

%     Version History:
%        Mar. 15, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check data structure
error(seizmocheck(data));
nrecs=numel(data);

% check phaselist
if(nargin==1); phaselist=[]; end
if(~isempty(phaselist) && ~isstring(phaselist) && ~iscellstr(phaselist))
    error('seizmo:arrivals2picks:badInput',...
        'PHASELIST must be a string like ''P'' or ''P,PP''!');
end

% default/check idx
if(nargin<3 || isempty(idx)); idx=0:9; end
if(~isnumeric(idx) || numel(idx)>10 ...
        || any(idx~=fix(idx) | idx<0 | idx>9))
    error('seizmo:arrivals2picks:badInput',...
        'IDX must be integers from 0-9 (eg. which KT/T/USER fields)!');
end
idx=unique(idx(:));

% toggle off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt insertion
try
    % get header info
    [o,kt,t,user]=getheader(data,'o','kt','t','user');
    
    % fix o time
    if(any(isnan(o)))
        warning('seizmo:arrivals2picks:unsetO',...
            'O header field unset for some records!  Treating as zero.');
    end
    o(isnan(o))=0;
    
    % get arrivals
    [arr,didx]=findarrivals(data,phaselist,varargin{:});
    
    % assign
    for i=1:nrecs
        if(any(i==didx))
            pick=find(i==didx);
            maxpick=min(numel(pick),numel(idx));
            t(i,idx(1:maxpick)+1)=[arr(pick(1:maxpick)).time]+o(i);
            kt(i,idx(1:maxpick)+1)={arr(pick(1:maxpick)).phase};
            user(i,idx(1:maxpick)+1)=[arr(pick(1:maxpick)).rayparameter];
        end
    end
    
    % update header
    data=changeheader(data,'t',t,'kt',kt,'user',user);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
