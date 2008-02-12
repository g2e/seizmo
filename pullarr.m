function [times]=pullarr(data,phase)
%PULLARR    Returns stored phase arrival time from SAClab data header
%
%    Description: Searches 'kt(n)' header fields for the specified phase. 
%     If found, returns the matching 't(n)' value.  Will return undefined
%     if a phase name match is not found.  In case of multiple entries,
%     will return only the first match found for each record - lower ktn/tn 
%     has preference.  NAME OF PHASE IS CASE SENSITIVE!
%
%    Usage: [times]=pullarr(data,phase);
%
%    Examples:
%     Ptimes=pullarr(data,'P');
%     sSKStimes=pullarr(data,'sSKS')
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: qcksnr, gh

% check nargin
error(nargchk(2,2,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% number of seismograms
nrecs=length(data);

% preallocate output
times=zeros(nrecs,1);
    
% grab header values
[kt,t]=gh(data,'kt','t');

% remove spaces from kt
kt=strtrim(kt);

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% do operations individually
for i=1:nrecs
    % find first match
    pos=find(strcmp(phase,kt(i,:)),1);
    
    % check for failure
    if(isempty(pos))
        warning('SAClab:noPhase','Could not phase %s for record %d',phase,i);
        v=data(i).version==vers;
        times(i)=h(v).undef.ntype;
        continue;
    end
    
    % add time
    times(i)=t(i,pos);
end

end