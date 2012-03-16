function [picks,idx]=findpicks(data,phaselist,varargin)
%FINDPICKS    Returns specified picks from SEIZMO records
%
%    Usage:    picks=findpicks(data,phaselist)
%              picks=findpicks(data,phaselist,k)
%              picks=findpicks(data,phaselist,k,'last')
%              [picks,idx]=findpicks(...)
%
%    Description:
%     PICKS=FINDPICKS(DATA,PHASELIST) returns the phase travel times for
%     the specified phases in PHASELIST stored in the pick header fields
%     (KT & T) of SEIZMO struct DATA.  PHASELIST may be either a cell array
%     of phases or a comma-separated string.  Note that phases are returned
%     in the order they are stored in the header, not in the order from
%     PHASELIST.
%
%     PICKS=FINDPICKS(DATA,PHASELIST,K) specifies the maximum number of
%     phases to return for each record.  Returns the first K picks
%     found.
%
%     PICKS=FINDPICKS(DATA,PHASELIST,K,'LAST') returns the last K picks
%     found matching the PHASELIST.
%
%     [PICKS,IDX]=FINDPICKS(...) also returns the indices indicating from
%     which field each pick in PICKS came.  This is particularly useful for
%     setting the IZTYPE header field (eg with TIMESHIFT).
%
%    Notes:
%
%    Examples:
%     % Only return the first P pick:
%     p=findpicks(data,'P',1);
%
%     % Shifting to the 1st P arrival and setting IZTYPE:
%     [P,i]=findpicks(data,'P',1);
%     data=timeshift(data,-P,strcat('it',num2str(i)));
%
%    See also: MAKEARRIVALS, ARRIVALS2PICKS, FINDARRIVALS

%     Version History:
%        Mar. 15, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check data structure
error(seizmocheck(data));
nrecs=numel(data);

% check phaselist
if(isstring(phaselist))
    phaselist=getwords(phaselist,',');
end
if(~iscellstr(phaselist))
    error('seizmo:findpicks:badInput',...
        'PHASELIST must be a string like ''P'' or ''P,PP''!');
end

% grab header values
[kt,t]=getheader(data,'kt','t');

% preallocation
picks=nan(nrecs,10);
idx=nan(nrecs,10);
n=zeros(nrecs,1);
maxpicks=1;

% find specified phases
for i=1:nrecs
    ok=find(ismember(kt(i,:),phaselist),varargin{:});
    n(i)=numel(ok);
    if(n(i)>0)
        if(n(i)>maxpicks); maxpicks=n(i); end
        picks(i,1:n(i))=t(i,ok);
        idx(i,1:n(i))=ok;
    end
end

% trim output
picks=picks(:,1:maxpicks);
idx=idx(:,1:maxpicks);

% throw warning for missed picks
if(any(n==0))
    if(all(n==0))
        tmp=['ALL OF THEM: 1 to ' num2str(nrecs)];
    elseif(sum(n==0)>20)
        tmp=[sprintf('%d ',find(n==0,10,'first')) ...
            '... ' sprintf('%d ',find(n==0,10,'last'))];
    else
        tmp=sprintf('%d ',find(n==0));
    end
    warning('seizmo:findpicks:noPhase',...
        ['Could not find phase(s):\n ' sprintf('%s ',phaselist{:}) ...
        '\nRecord(s):\n ' tmp]);
end

end
