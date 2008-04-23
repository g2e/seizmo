function [data]=rms(data,nsamples)
%RMS    Return sliding window root-mean-square of SAClab data records
%
%    Description: DATA_OUT=RMS(DATA_IN,NSAMPLES) takes the sliding window
%     root-mean-square of the records in DATA_IN using a window of length 
%     NSAMPLES samples centered on each point in the original records.  
%     Windows extending outside the record are padded with zeros, which 
%     gives the output tapered edges. NSAMPLES can be a scalar (each
%     record has the same window size) or a vector (define each record's
%     window size separately).  Even NSAMPLES will be changed to NSAMPLES+1
%     so that the sliding window is centered on an existing time point.
%
%    Usage: [data]=rms(data,nsamples)
%
%    Examples:
%     
%      Compare envelope and a 10-sample sliding window rms:
%       p2([envelope(data(1)) rms(data(1),10)])
%
%    See also: envelope, robustrms

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || ~isvector(nsamples) || any(nsamples<1) ...
        || length(nsamples)~=nrecs || any(fix(nsamples)~=nsamples))
    error('SAClab:rms:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% half window length
half1=floor(nsamples/2);

% sliding window rms
for i=1:nrecs
    oclass=str2func(class(data(i).x));
    
    % filter to avoid loop
    data(i).x=oclass(sqrt(filter(ones(1,nsamples(i))/nsamples(i),1,...
        [double(data(i).x).^2; zeros(half1(i),size(data(i).x,2))])));
    
    % centered windows
    data(i).x(1:half1(i),:)=[];
    
    % update header
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
