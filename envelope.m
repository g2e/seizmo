function [data]=envelope(data,n)
%ENVELOPE    Return envelope of SAClab data records
%
%    Description: Calculates and returns the envelopes of the data records.
%      The envelope is the complex magnitude of a records analytic signal.
%     Uses the Matlab function hilbert (Signal Processing Toolbox).
%
%    Usage: [data]=envelope(data);
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: hilbrt

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
        % get envelope
        oclass=str2func(class(data(i).x));
        data(i).x=abs(hilbert(double(data(i).x),n));
        data(i).x=oclass(data(i).x);
        
        % update header
        data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
            'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
    end
end

end