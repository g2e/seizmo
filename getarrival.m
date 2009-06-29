function [times]=getarrival(data,phase)
%GETARRIVAL    Returns stored phase arrival time from SEIZMO data header
%
%    Usage:    times=getarrival(data,phase)
%
%    Description: GETARRIVAL(DATA,PHASE) searches 'kt(n)' header fields in
%     the SEIZMO structure DATA for the specified phase PHASE.  If found,
%     the matching 't(n)' value is returned.  If not, NaN is returned. In
%     case of multiple entries for the same phase, only the first match
%     found for each record is returned - lower index has preference.  Note
%     that the returned time is based on the reference time and is not
%     relative to the origin (see example below).
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
%        June 29, 2009 - add testing table, doc update
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
%     Last Updated June 29, 2009 at 04:30 GMT

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
times=nan(nrecs,1);
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
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
