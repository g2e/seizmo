function [data]=rms(data,nsamples)
%RMS    Return sliding-window root-mean-square of SAClab data records
%
%    Description: RMS(DATA,NSAMPLES) takes the sliding window
%     root-mean-square of the records in DATA using a window of length 
%     NSAMPLES samples centered on each point in the original records.  
%     Windows extending outside the record are padded with zeros, which 
%     gives the output tapered edges. NSAMPLES can be a scalar (each
%     record has the same window size) or a vector (define each record's
%     window size separately).  
%
%    Notes:
%     - If NSAMPLES is even it will be changed to NSAMPLES+1 so that the 
%       sliding window is centered on an existing time point.
%
%    System requirements: Matlab 7
%
%    Data requirements: None
%
%    Usage: data=rms(data,nsamples)
%
%    Examples:
%      Compare an envelope and a 10-sample sliding window rms:
%       p2([envelope(data(1)) rms(data(1),10)])
%
%    See also: envelope, robustrms

%     Version History:
%        Apr.  9, 2008 - initial version
%        Apr. 23, 2008 - changed behavior for windows with even npts
%        May  12, 2998 - dep* fix
%        July 17, 2008 - history update, documentation update, .dep rather
%                        than .x, dataless handling
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 18, 2008 at 00:55 GMT

% todo:
%

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'dep'))

% number of records
nrecs=numel(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || any(nsamples<1) || numel(nsamples)~=nrecs ...
        || any(fix(nsamples)~=nsamples))
    error('SAClab:rms:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% half window length
half1=floor(nsamples/2);

% sliding window rms
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    
    % filter to avoid loop
    data(i).dep=oclass(sqrt(filter(ones(1,nsamples(i))/nsamples(i),1,...
        [double(data(i).dep).^2; zeros(half1(i),size(data(i).dep,2))])));
    
    % center windows
    data(i).dep(1:half1(i),:)=[];
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
