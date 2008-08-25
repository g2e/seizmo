function [data]=robustrms(data,nsamples)
%ROBUSTRMS    Return sliding-window root-median-square of SAClab data records
%
%    Description: ROBUSTRMS(DATA,NSAMPLES) takes the sliding window
%     root-median-square of the records in DATA using a window of length 
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
%    Usage: data=robustrms(data,nsamples)
%
%    Examples:
%      A comparison of envelope, rms and robustrms:
%       p2([envelope(data(1)) rms(data(1),10) robustrms(data(1),10)])
%
%    See also: envelope, rms

%     Version History:
%        Apr.  9, 2008 - initial version
%        Apr. 23, 2008 - changed behavior for windows with even npts
%        May  12, 2998 - dep* fix
%        July 18, 2008 - history update, documentation update, .dep rather
%                        than .x, dataless handling
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 18, 2008 at 19:45 GMT

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
    error('SAClab:robustrms:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% half window length
half1=floor(nsamples/2);
half12=2*half1;

% sliding window rms
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    
    % pad record with zeros and square it
    sz=size(data(i).dep);
    data(i).dep=[zeros(half1(i),sz(2)); ...
                double(data(i).dep).^2; ...
                zeros(half1(i),sz(2))];
    
    % get indices of windows
    indices=hankel(1:sz(1),sz(1):sz(1)+half12(i)).';
    indices=repmat(indices,[1 1 sz(2)])...
        +repmat(permute((0:sz(2)-1)*(sz(1)+half12(i)),[1 3 2]),...
            [nsamples(i) sz(1) 1]);
    
    % get root-median
    data(i).dep=...
        oclass(sqrt(permute(median(data(i).dep(indices),1),[2 3 1])));
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
