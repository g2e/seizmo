function [arr,idx]=findarrivals(data,phaselist,varargin)
%FINDARRIVALS    Returns specified arrival info from SEIZMO records
%
%    Usage:    arr=findarrivals(data)
%              arr=findarrivals(data,phaselist)
%              arr=findarrivals(data,phaselist,k)
%              arr=findarrivals(data,phaselist,k,'last')
%              [arr,idx]=findarrivals(...)
%
%    Description:
%     ARR=FINDARRIVALS(DATA) returns the arrival info stored in SEIZMO
%     struct DATA.  See MAKEARRIVALS & TAUPTIME for more info.  ARR is a
%     struct formatted the same as TAUPTIME output.
%
%     ARR=FINDARRIVALS(DATA,PHASELIST) returns the arrival info for the
%     specified phases in PHASELIST.  PHASELIST may be either a cell array
%     of phases or a comma-separated string.  If PHASELIST is empty ([])
%     then all arrivals are returned.
%
%     ARR=FINDARRIVALS(DATA,PHASELIST,K) specifies the maximum number of
%     phases to return for each record.  Returns the first K arrivals
%     found.
%
%     ARR=FINDARRIVALS(DATA,PHASELIST,K,'LAST') returns the last K arrivals
%     found matching the PHASELIST.
%
%     [ARR,IDX]=FINDARRIVALS(...) also returns the indices indicating from
%     which record each arrival in ARR came.
%
%    Notes:
%     - Calls MAKEARRIVALS if some records are missing .misc.arrivals
%
%    Examples:
%     % Only return the first P arrival:
%     arr=findarrivals(data,'P',1);
%
%    See also: MAKEARRIVALS, ARRIVALS2PICKS, FINDPICKS

%     Version History:
%        Mar. 15, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% check data structure
error(seizmocheck(data));
nrecs=numel(data);

% check phaselist
if(nargin==1); phaselist=[]; end
if(~isempty(phaselist) && ~isstring(phaselist) && ~iscellstr(phaselist))
    error('seizmo:findarrivals:badInput',...
        'PHASELIST must be a string like ''P'' or ''P,PP''!');
end

% find records that need makearrivals called
call=false(nrecs,1);
for i=1:nrecs
    if(~isfield(data(i).misc,'arrivals'))
        call(i)=true;
    end
end

% call makearrivals if needed
if(any(call))
    try
        oldseizmocheckstate=seizmocheck_state(false);
        data(call)=makearrivals(data(call),phaselist);
        seizmocheck_state(oldseizmocheckstate);
    catch
        seizmocheck_state(oldseizmocheckstate);
        error(lasterror);
    end
end

% parse comma-separated phases
if(isstring(phaselist))
    phaselist=getwords(phaselist,',');
end

% support for return all via empty phaselist
if(isempty(phaselist))
    if(nargin>2 && isscalar(varargin{1}) && isnumeric(varargin{1}) ...
        && varargin{1}==fix(varargin{1}) && varargin{1}>0)
        k=varargin{1};
    else
        k=[];
    end
else
    k=[];
end

% find specified phases
arr=[]; idx=[];
none=false(nrecs,1);
for i=1:nrecs
    phases={data(i).misc.arrivals.phase}';
    if(isempty(phaselist))
        if(isempty(k))
            ok=1:numel(phases);
        else
            ok=1:min(numel(phases),k);
        end
    else
        ok=find(ismember(phases,phaselist),varargin{:});
    end
    none(i)=isempty(ok);
    if(~none(i))
        arr=cat(1,arr,data(i).misc.arrivals(ok));
        idx=cat(1,idx,i*ones(numel(ok),1));
    end
end

% throw warning for missed arrivals
if(any(none))
    if(all(none))
        tmp=['ALL OF THEM: 1 to ' num2str(nrecs)];
    elseif(sum(none)>20)
        tmp=[sprintf('%d ',find(none,10,'first')) ...
            '... ' sprintf('%d ',find(none,10,'last'))];
    else
        tmp=sprintf('%d ',find(none));
    end
    warning('seizmo:findarrivals:noPhase',...
        ['Could not find phase(s):\n ' sprintf('%s ',phaselist{:}) ...
        '\nRecord(s):\n ' tmp]);
end

end
