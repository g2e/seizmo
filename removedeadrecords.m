function [data]=removedeadrecords(data,option)
%REMOVEDEADRECORDS    Removes constant SEIZMO records
%
%    Description: REMOVEDEADRECORDS(DATA) removes records that have no
%     change in the dependent component.  These can cause problems in
%     analysis and are not worth keeping.  Uses the header fields 
%     depmin/depmax so that records can be eliminated before reading
%     in the data.
%
%     REMOVEDEADRECORDS(DATA,OPTION) allows changing how records are
%     determined as dead.  OPTION is a logical that when set true
%     (default), will use the header fields depmin/depmax to look for dead
%     records.  When OPTION is false, the data (in .dep) is used to find
%     dead records.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    data=removedeadrecords(data)
%
%    Examples:
%     Remove dead records before reading in data from current directory:
%      data=readdata(removedeadrecords(readheaders('*')));
%
%    See also: removemean, removetrend

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 28, 2008 - SEISCHK support
%        Mar.  4, 2008 - minor doc update
%        Nov. 22, 2008 - update for new name schema (now REMOVEDEADRECORDS)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 08:00 GMT

% todo:

% check input
error(nargchk(1,2,nargin))

% check option
if(nargin==1 || isempty(option))
    option=true;
elseif(~islogical(option) || ~isscalar(option))
    error('seizmo:removedeadrecords','OPTION must be true or false!');
end

% proceed by option
if(option)
    % check data structure
    error(seizmocheck(data))
    
    % turn off struct checking
    oldseizmocheckstate=get_seizmocheck_state;
    set_seizmocheck_state(false);

    % check headers
    data=checkheader(data);
    
    % get depmax, depmin
    [depmax,depmin]=getheader(data,'depmax','depmin');
    
    % remove dead records
    data(depmax-depmin==0)=[];
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
else
    % check data structure
    error(seizmocheck(data,'dep'))
    
    % get min/max of data
    nrecs=numel(data);
    dmin=nan(nrecs,1); dmax=dmin;
    for i=1:nrecs
        dmin(i)=min(data(i).dep(:));
        dmax(i)=max(data(i).dep(:));
    end
    
    % remove dead records
    data((dmax-dmin)==0)=[];
end

end
