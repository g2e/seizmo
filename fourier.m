function [data]=fourier(data,format,pp2)
%FOURIER    Converts time domain SAClab records to the frequency domain
%
%    Description: Converts SAClab records from the time domain to the
%     frequency domain using the fast fourier transform.  Following SAC
%     formatting, an output choice between real-imaginary and amplitude-
%     phase is allowed through the 'format' input ('amph' or 'rlim' only -
%     default is 'amph').  The power of 2 padding can be adjusted with the
%     pp2 input (fft length = 2^(nextpow2(len)+pp2) , default pp2=1).
%
%    Note:
%     - SAC (and thus SAClab for sanity) stores amp/real/imag spectral 
%       data in a raw form that first has to be scaled to give correct 
%       amplitudes.  Divide records (except for phase!) by npts*delta/2 to 
%       get accurate spectral information.  All SAClab functions that work
%       with fourier will expect that this scaling has NOT been applied.
%
%    Usage: data=fourier(data,format,pp2)
%
%    Examples:
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
    error('SAClab:fourier:badInput','bad format string: %s',format)
end
if(~isnumeric(pp2) || ~isscalar(pp2) || fix(pp2)~=pp2)
    error('SAClab:fourier:badInput','pp2 must be a scalar integer')
end

% retreive header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
[b,delta]=gh(data,'b','delta');

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('SAClab:fourier:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Time Series File'))...
        && any(~strcmp(iftype,'General X vs Y file')))
    error('SAClab:fourier:illegalOperation',...
        'illegal filetype for this operation')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of datum/components
    [len,ncmp]=size(data(i).x);
    
    % get frequency info
    nspts=2^(nextpow2(len)+pp2);
    sb=0; se=1/(delta(i)*2); sdelta=2*se/nspts;
    
    % fft
    data(i).x=delta(i)*fft(data(i).x,nspts);    % SAC compatible
    %data(i).x=2*fft(data(i).x,nspts)/len;      % True amplitudes
    
    % expand data to make room for split
    data(i).x(:,2:2:2*ncmp)=data(i).x;
    
    % split complex by desired filetype
    if(strcmpi(format,'rlim'))
        data(i).x(:,1:2:end)=real(data(i).x(:,2:2:end));
        data(i).x(:,2:2:end)=imag(data(i).x(:,2:2:end));
        data(i)=ch(data(i),'iftype','Spectral File-Real/Imag');
    else
        data(i).x(:,1:2:end)=abs(data(i).x(:,2:2:end));
        data(i).x(:,2:2:end)=angle(data(i).x(:,2:2:end));
        data(i)=ch(data(i),'iftype','Spectral File-Ampl/Phase');
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header (note there is no field 'se')
    data(i)=ch(data(i),'b',sb,'e',se,'delta',sdelta,'sb',b(i),...
        'sdelta',delta(i),'nspts',len,'npts',nspts);
end

end
