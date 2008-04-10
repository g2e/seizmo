function [data]=robustrms(data,nsamples)
%ROBUSTRMS    Return sliding window root-median-square of SAClab data records
%
%    Description: DATA_OUT=RMS(DATA_IN,NSAMPLES) takes the sliding window
%     root-median-square of the records in DATA_IN using a window of length 
%     NSAMPLES samples centered on each point in the original records.  
%     Windows extending outside the record are padded with zeros, which 
%     gives the output tapered edges. NSAMPLES can be a scalar (each
%     record has the same window size) or a vector (define each record's
%     window size separately).
%
%    Usage: [data]=rms(data,nsamples)
%
%    Examples:
%     
%      A comparison of envelope, rms and robust rms:
%       p2([envelope(data(1)) rms(data(1),10) robustrms(data(1),10)])
%
%    See also: envelope, rms

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% number of records
nrecs=length(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || ~isvector(nsamples) ...
        || length(nsamples)~=nrecs || any(fix(nsamples)~=nsamples) ...
        || any(nsamples<1) || any(~mod(nsamples,2)))
    error('SAClab:rms:badInput',...
        'NSAMPLES must be a scalar/vector of odd positive integer(s)')
end

% half window length
half1=floor(nsamples/2);
half12=2*half1;

% sliding window rms
for i=1:nrecs
    % get storage class of data
    oclass=str2func(class(data(i).x));
    
    % pad record with zeros and square it
    sz=size(data(i).x);
    data(i).x=[zeros(half1(i),sz(2)); ...
                double(data(i).x).^2; ...
                zeros(half1(i),sz(2))];
    
    % get indices of windows
    indices=hankel(1:sz(1),sz(1):sz(1)+half12(i)).';
    indices=repmat(indices,[1 1 sz(2)])...
        +repmat(permute((0:sz(2)-1)*(sz(1)+half12(i)),[1 3 2]),...
            [nsamples(i) sz(1) 1]);
    
    % get root-median
    data(i).x=oclass(sqrt(permute(median(data(i).x(indices),1),[2 3 1])));
    
    % update header
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
