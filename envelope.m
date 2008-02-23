function [data]=envelope(data,n)
%ENVELOPE    Return envelope of SAClab data records
%
%    Description: Returns the envelopes of the data records, which is the
%     complex magnitude of a record's analytic signal.  Uses the Matlab 
%     function hilbert (Signal Processing Toolbox).
%
%    Usage: [data]=envelope(data);
%
%    See also: hilbrt

% check nargin
error(nargchk(1,2,nargin))

% default fft length to none (lets fft do it)
if(nargin==1); n=[]; end

% check data structure
if(~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% check spacing
if(any(~strcmp(glgc(data,'leven'),'true')))
    error('SAClab:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced data');
end

% do operations individually
for i=1:length(data)
    oclass=str2func(class(data(i).x));
    data(i).x=abs(hilbert(double(data(i).x),n));
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end

