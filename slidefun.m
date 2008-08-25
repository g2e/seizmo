function [data]=slidefun(data,fun,nsamples)
%SLIDEFUN    Apply a sliding window function to SAClab data records
%
%    Description: SLIDEFUN(DATA,FUN,NSAMPLES) applies a function defined by
%     the function handle FUN to a sliding window of NSAMPLES samples to 
%     the dependent component(s) of SAClab data records in DATA.
%     
%     FUN is expected to handle a column vector (single component) or an 
%     array (multiple component - components are distributed in columns).  
%     Output of FUN is assigned to the time point corresponding to the 
%     center of the sliding window.  Vector output will be distributed as 
%     multiple components corresponding to a single time.  Array output 
%     will produce an error.
%     
%     NSAMPLES can be a scalar (each record has the same window size) or a 
%     vector (define each record's window size separately).  If NSAMPLES is
%     even, it will be changed to NSAMPLES+1 so that the window will be
%     centered on an existing time point.  Records are zero padded so that
%     the first window is centered on the first point in the record and the 
%     last window is centered on the last point in the record.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record.
%
%    System requirements: Matlab 7
%
%    Data requirements: NONE
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=slidefun(data,fun,nsamples)
%
%    Examples:
%     Running absolute mean from G. D. Bensen et al, 2006:
%      slidefun(data,@(x)mean(abs(x)),ceil(1/(2*delta*fmin)))
%
%    See also: seisfun, rms, robustrms

%     Version History:
%        Apr. 23, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 18, 2008 - documentation update, dataless support, .dep
%                        rather than .x, added history, single ch call
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 18, 2008 at 23:15 GMT

% todo:
%

% check nargin
error(nargchk(3,3,nargin))

% check data structure
error(seischk(data,'dep'))

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('SAClab:slidefun:badInput','FUN must be a function handle!')
end

% number of records
nrecs=numel(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || any(nsamples<1) || numel(nsamples)~=nrecs ...
        || any(fix(nsamples)~=nsamples))
    error('SAClab:slidefun:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% half window length
half1=floor(nsamples/2);
half12=2*half1;

% loop through each record
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    
    % pad record with zeros
    sz=size(data(i).dep);
    data(i).dep=[zeros(half1(i),sz(2)); ...
                double(data(i).dep); ...
                zeros(half1(i),sz(2))];
    
    % sliding window for each point
    % (starts at last point to allocate temp once)
    for j=sz(1):-1:1
        temp(j,:)=fun(data(i).dep(j:j+half12,:)); %#ok<AGROW>
    end
    
    % assign back to data and change storage back
    data(i).dep=oclass(temp);
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
