function [data]=envelope(data,pp2)
%ENVELOPE    Return envelopes of SAClab data records
%
%    Description: Returns the envelopes of the SAClab data records (complex
%     magnitude of a record's analytic signal).  Uses the Matlab function
%     hilbert (Signal Processing Toolbox).
%
%    Usage: [data]=envelope(data)
%
%    Examples:
%
%    See also: hilbrt

% check nargin
error(nargchk(1,2,nargin))

% check data structure
error(seischk(data,'x'))

% default fft length to next power of 2 + 1
if(nargin==1); pp2=1; end

% check spacing
if(any(~strcmp(glgc(data,'leven'),'true')))
    error('SAClab:envelope:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced data');
end

% get envelopes
for i=1:length(data)
    len=size(data(i).x,1);
    nfft=2^(nextpow2(len)+pp2);
    oclass=str2func(class(data(i).x));
    data(i).x=abs(hilbert(double(data(i).x),nfft));
    data(i).x=oclass(data(i).x(1:len,:));
    
    % update header
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
