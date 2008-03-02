function [data]=fourier(data,format,pp2)
%FOURIER    Converts time domain seislab records to the frequency domain
%
%    Description: Converts seislab records from the time domain to the
%     frequency domain using the fast fourier transform.  Following SAC
%     formatting, an output choice between real-imaginary and amplitude-
%     phase is allowed through the 'format' input ('amph' or 'rlim' only -
%     default is 'amph').  The power of 2 padding can be adjusted with the
%     pp2 input (fft length = 2^(nextpow2(len)+pp2) , default pp2=1).
%
%    Note:
%     - Amplitudes in the 'amph' case will not match those from SAC.  The
%       convention here is 1/npts times the SAC amplitudes.  This
%       convention gives actual spectral amplitudes.
%
%    Usage: data=fft(data,format,pp2)
%
%    See also: ifourier

% check nargin
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'x'))

% defaults
if(nargin<3 || isempty(pp2)); pp2=1; end
if(nargin<2 || isempty(format)); format='amph'; end

% check inputs
if(~any(strcmpi(format,{'amph' 'rlim'})))
    error('seislab:fourier:badInput','bad format string: %s',format)
end
if(~isnumeric(pp2) || ~isscalar(pp2) || fix(pp2)~=pp2)
    error('seislab:fourier:badInput','pp2 must be a scalar integer')
end

% retreive header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
[b,delta]=gh(data,'b','delta');

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('seislab:fourier:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Time Series File'))...
        && any(~strcmp(iftype,'General X vs Y file')))
    error('seislab:fourier:illegalOperation',...
        'illegal filetype for this operation')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of datum/components
    ncmp=size(data(i).x);
    
    % get frequency info
    nspts=2^(nextpow2(ncmp(1))+pp2);
    sb=0; se=1/(delta(i)*2); sdelta=2*se/nspts;
    
    % fft
    data(i).x=fft(data(i).x,nspts)/ncmp(1);
    
    % split complex by desired filetype
    if(strcmpi(format,'rlim'))
        data(i).x(:,2)=imag(data(i).x(:,1));
        data(i).x(:,1)=real(data(i).x(:,1));
        data(i)=ch(data(i),'iftype','Spectral File-Real/Imag');
    else
        data(i).x(:,2)=angle(data(i).x(:,1));
        data(i).x(:,1)=2*abs(data(i).x(:,1));
        data(i)=ch(data(i),'iftype','Spectral File-Ampl/Phase');
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header (note there is no field 'se')
    data(i)=ch(data(i),'b',sb,'e',se,'delta',sdelta,'sb',b(i),...
        'sdelta',delta(i),'nspts',ncmp(1),'npts',nspts);
end

end
