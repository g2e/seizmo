function [times,n]=getarrival(data,phase)
%GETARRIVAL    Returns stored phase arrival time from SEIZMO data header
%
%    Usage:    times=getarrival(data,phase)
%              [times,n]=getarrival(data,phase)
%
%    Description: GETARRIVAL(DATA,PHASE) searches 'kt(n)' header fields in
%     the SEIZMO structure DATA for the specified phase PHASE.  If found,
%     the matching 't(n)' value is returned.  If not, NaN is returned. In
%     case of multiple entries for the same phase, only the first match
%     found for each record is returned - lower index has preference.  Note
%     that the returned time is based on the reference time and is not
%     relative to the origin (see example below).
%
%     [TIMES,N]=GETARRIVAL(DATA,PHASE) also returns an index array that
%     indicates the header field from which the time was found for each
%     record.  So [3 2 1] would mean record 1's time came from t3, record
%     2's time from t2 and record 3's from t1.  This is particularly useful
%     for setting iztype when combining getarrival with timeshift.
%     
%    Notes:
%     - NAME OF PHASE IS CASE SENSITIVE!
%     - Returned times are relative to the reference time!
%
%    Examples:
%     Ptimes=getarrival(data,'P');
%     sSKStimes=getarrival(data,'sSKS')
%
%     Get arrival time that is relative to origin time:
%      Ptimes=getarrival(data,'P')-getheader(data,'o');
%
%    See also: quicksnr, getheader, addarrival

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 29, 2008 - uses SEISCHK now
%        Mar.  4, 2008 - minor doc update
%        Oct.  8, 2008 - doc update, add history
%        Nov. 16, 2008 - doc update, history fix, rename from PULLARR to
%                        GETARRIVAL, minor code clean
%        Nov. 24, 2008 - minor code cleaning
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 29, 2009 - doc update
%        June 30, 2009 - second output: t index
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:55 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% toggle off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% grab header values
[kt,t]=getheader(data,'kt','t');

% remove spaces from kt
kt=strtrim(kt);

% number of records
nrecs=numel(data);

% do operations individually
times=nan(nrecs,1); n=times;
for i=1:nrecs
    % find first match
    pos=find(strcmp(phase,kt(i,:)),1);
    
    % check for failure
    if(isempty(pos))
        warning('seizmo:getarrival:noPhase',...
            'Could not find phase %s for record %d !',phase,i);
        continue;
    end
    
    % add time
    times(i)=t(i,pos);
    n(i)=pos-1;
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
