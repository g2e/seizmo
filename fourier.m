function [data]=fourier(data,format,pp2)
%FOURIER    Converts timeseries SAClab records to the frequency domain
%
%    Description: Converts SAClab records from the time domain to the
%     frequency domain using the fast fourier transform.  Following SAC
%     formatting, an output choice between real-imaginary and amplitude-
%     phase is allowed through the 'format' input ('amph' or 'rlim' only -
%     default is 'amph').  The power of 2 padding can be adjusted with the
%     pp2 input (fft length = 2^(nextpow2(len)+pp2) , default pp2=1).
%
%    Usage: data=fft(data,format,pp2)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: ifourier, selctcmp

% check nargin
error(nargchk(1,3,nargin))

% defaults
if(nargin<3 || isempty(pp2)); pp2=1; end
if(nargin<2 || isempty(format)); format='amph'; end

% check
if(~any(strcmpi(format,{'amph' 'rlim'})))
    error('bad output format string: %s',format)
end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% retreive header info
[b,delta,leven,iftype]=gh(data,'b','delta','leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% loop through records
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of components
    ncmp=size(data(i).x);
    
    % check spacing and filetype
    if(leven(i)~=h(v).true)
        warning('SAClab:illegalOperation','Illegal operation on unevely spaced record %d',i);
        continue;
    elseif(iftype(i)==h(v).enum(1).val.ixyz)
        warning('SAClab:illegalOperation','Illegal operation on xyz file: record %d',i);
        continue;
    elseif(iftype(i)==h(v).enum(1).val.irlim || iftype(i)==h(v).enum(1).val.iamph)
        warning('SAClab:illegalOperation','Illegal operation on spectral file: record %d',i);
        continue;
    elseif(isfield(h(v).enum(1).val,'incmp') && iftype(i)==h(v).enum(1).val.incmp)
       warning('SAClab:illegalOperation','Illegal operation on multi-component file: record %d',i);
        continue;
    elseif(ncmp(2)>1)
       warning( 'SAClab:illegalOperation','Illegal operation on multi-component record %d',i);
       continue;
    end
    
    % get frequency info
    nspts=2^(nextpow2(ncmp(1))+pp2);
    sb=0; se=1/(delta(i)*2); sdelta=2*se/nspts;
    
    % fft
    data(i).x=fft(data(i).x,nspts)/ncmp(1);
    
    % split complex by desired filetype
    if(strcmpi(format,'rlim'))
        data(i).x(:,2)=imag(data(i).x);
        data(i).x(:,1)=real(data(i).x(:,1));
        data(i)=ch(data(i),'iftype','irlim');
    else
        data(i).x(:,2)=angle(data(i).x);
        data(i).x(:,1)=2*abs(data(i).x(:,1));
        data(i)=ch(data(i),'iftype','iamph');
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header (note there is no field 'se')
    data(i)=ch(data(i),'b',sb,'e',se,'delta',sdelta,'sb',b(i),...
        'sdelta',delta(i),'nspts',ncmp(1),'npts',nspts);
end

end