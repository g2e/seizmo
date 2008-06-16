function [data]=hilbrt(data,pp2)
%HILBRT    Return Hilbert transform of SAClab data records
%
%    Description: Calculates and returns the Hilbert transform of SAClab
%     data records.  Hilbert tranform adds a 90 degree phase shift.  Uses
%     the Matlab function hilbert (Signal Processing Toolbox).
%
%    Usage: [data]=hilbrt(data)
%
%    Examples:
%     To do a negative 90 degree phase shift:
%      data=hilbert(mul(data,-1))
%
%    See also: envelope, mul

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

% get Hilbert transforms
for i=1:length(data)
    len=size(data(i).x,1);
    nfft=2^(nextpow2(len)+pp2);
    oclass=str2func(class(data(i).x));
    data(i).x=imag(hilbert(double(data(i).x),nfft));
    data(i).x=oclass(data(i).x(1:len,:));
    
    % update header
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
