function [data]=slidefun(data,fun,nsamples)
%SLIDEFUN    Apply a sliding window function to SAClab data records
%
%    Description: SLIDEFUN(DATA,FUN,NSAMPLES) applies a function FUN to a
%     sliding window of NSAMPLES samples to SAClab data records in DATA.  
%     Function FUN will only modify the amplitudes of the records and not 
%     the length or timing.  The number of components in a multi-component
%     record can be changed.
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
%    Usage: [data]=slidefun(data,fun,nsamples)
%
%    Note: Will update header values depmen, depmin, depmax.
%
%    Examples:
%
%     Benson et al 2006:
%      slidefun(data,@(x)mean(abs(x)),fix(maxper/(2*delta)))
%
%    See also: seisfun, rms, robustrms

% check nargin
error(nargchk(3,3,nargin))

% check data structure
error(seischk(data,'x'))

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('FUN must be a function handle!')
end

% number of records
nrecs=length(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || ~isvector(nsamples) || any(nsamples<1) ...
        || length(nsamples)~=nrecs || any(fix(nsamples)~=nsamples))
    error('SAClab:slidefun:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% half window length
half1=floor(nsamples/2);
half12=2*half1;

% sliding window rms
for i=1:nrecs
    % get storage class of data
    oclass=str2func(class(data(i).x));
    
    % pad record with zeros
    sz=size(data(i).x);
    data(i).x=[zeros(half1(i),sz(2)); ...
                double(data(i).x); ...
                zeros(half1(i),sz(2))];
    
    % sliding window for each point
    for j=sz(1):-1:1
        temp(j,:)=fun(data(i).x(j:j+half12,:)); %#ok<AGROW>
    end
    
    % assign back to data and change storage back
    data(i).x=oclass(temp);
    
    % update header
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
