function [times]=pullarr(data,phase)
%PULLARR    Returns stored phase arrival time from SAClab data header
%
%    Description: Searches 'kt(n)' header fields for the specified phase. 
%     If found, returns the matching 't(n)' value.  Will return undefined
%     if a phase name match is not found.  In case of multiple entries,
%     will return only the first match found for each record - lower ktn/tn 
%     has preference.  NAME OF PHASE IS CASE SENSITIVE!
%
%    Usage: [times]=pullarr(data,phase)
%
%    Examples:
%     Ptimes=pullarr(data,'P');
%     sSKStimes=pullarr(data,'sSKS')
%
%    See also: qcksnr, gh

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data))

% number of records
nrecs=length(data);

% preallocate output
times=nan(nrecs,1);
    
% grab header values
[kt,t]=gh(data,'kt','t');

% remove spaces from kt
kt=strtrim(kt);

% do operations individually
for i=1:nrecs
    % find first match
    pos=find(strcmp(phase,kt(i,:)),1);
    
    % check for failure
    if(isempty(pos))
        warning('SAClab:pullarr:noPhase',...
            'Could not find phase %s for record %d',phase,i);
        continue;
    end
    
    % add time
    times(i)=t(i,pos);
end

end
