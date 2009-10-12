function [data,removed]=removedeadrecords(data,option)
%REMOVEDEADRECORDS    Removes constant SEIZMO records
%
%    Usage:    data=removedeadrecords(data)
%              data=removedeadrecords(data,option)
%              [data,removed]=removedeadrecords(...)
%
%    Description: REMOVEDEADRECORDS(DATA) removes records that have no
%     change in the dependent component.  These can cause problems in
%     analysis and are not worth keeping.  Uses the header fields 
%     'depmin'/'depmax' so that records can be eliminated before reading
%     in the data.
%
%     REMOVEDEADRECORDS(DATA,OPTION) allows changing how records are
%     determined as dead.  OPTION is a logical that when set true
%     (default), will use the header fields depmin/depmax to look for dead
%     records.  When OPTION is false, the data (in .dep) is used to find
%     dead records.
%
%     [DATA,REMOVED]=REMOVEDEADRECORDS(...) also returns a listing of the
%     indices of the records removed in REMOVED.  These indices are
%     relative to the input dataset, not the output dataset (obviously
%     because the records are no longer in the output dataset!).
%
%    Notes:
%
%    Header changes: NONE (may use CHECKHEADER)
%
%    Examples:
%     Remove dead records before reading in data from current directory:
%      data=readdata(removedeadrecords(readheaders('*')));
%
%    See also: REMOVEMEAN, REMOVETREND

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 28, 2008 - SEISCHK support
%        Mar.  4, 2008 - minor doc update
%        Nov. 22, 2008 - update for new name schema (now REMOVEDEADRECORDS)
%        Dec.  8, 2008 - 2nd output lists removed indices
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check input
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check option
if(nargin==1 || isempty(option))
    option=true;
elseif(~islogical(option) || ~isscalar(option))
    error('seizmo:removedeadrecords','OPTION must be true or false!');
end

% proceed by option
if(option)
    % check data structure
	msg=seizmocheck(data);
	if(~isempty(msg)); error(msg.identifier,msg.message); end
    
    % turn off struct checking
    oldseizmocheckstate=get_seizmocheck_state;
    set_seizmocheck_state(false);

    % check headers
    data=checkheader(data);
    
    % get depmax, depmin
    [depmax,depmin]=getheader(data,'depmax','depmin');
    
    % remove dead records
    removed=((depmax-depmin)==0);
    data(removed)=[];
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
else
    % check data structure
	msg=seizmocheck(data,'dep');
	if(~isempty(msg)); error(msg.identifier,msg.message); end
    
    % get min/max of data
    nrecs=numel(data);
    dmin=nan(nrecs,1); dmax=dmin;
    for i=1:nrecs
        dmin(i)=min(data(i).dep(:));
        dmax(i)=max(data(i).dep(:));
    end
    
    % remove dead records
    removed=((dmax-dmin)==0);
    data(removed)=[];
end

end
