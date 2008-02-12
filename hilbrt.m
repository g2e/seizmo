function [data]=hilbrt(data,n)
%HILBRT    Return Hilbert transform of SAClab data records
%
%    Description: Calculates and returns the Hilbert transform of the data
%     records.  Hilbert tranform adds a 90 degree phase shift to the data.
%     Uses the Matlab function hilbert (Signal Processing Toolbox).
%
%    Usage: [data]=hilbrt(data);
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: envelope

% check nargin
error(nargchk(1,2,nargin))

% default fft length to none (let fft do it)
if(nargin==1); n=[]; end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% header info
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end
leven=gh(data,'leven');

% do operations individually
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % must be evenly spaced
    if(leven(i)~=h(v).true)
        error('SAClab:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    else
        % get Hilbert transformation
        oclass=str2func(class(data(i).x));
        data(i).x=imag(hilbert(data(i).x,n));
        data(i).x=oclass(data(i).x);
        
        % update header
        data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
            'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
    end
end

end