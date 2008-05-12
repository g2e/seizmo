function [data]=dft(data,format,pow2pad)
%DFT    Performs a discrete fourier transform on SAClab data records
%
%    Description: DFT(DATA,FORMAT) converts SAClab records from the time 
%     domain to the frequency domain using a discrete fourier transform
%     (Matlab's fft).  Following SAC formatting, an option FORMAT can be
%     used to select between real-imaginary and amplitude-phase formats
%     (FORMAT can be either 'amph' or 'rlim' - default is 'amph').  
%
%     DFT(DATA,FORMAT,POW2PAD) lets the power of 2 zero padding be adjusted
%     using POW2PAD (default value is 1) according to the formula:
%       fftlength=2^(nextpow2(NPTS)+POW2PAD)
%
%     In the frequency domain B, DELTA, and NPTS are changed to the 
%     beginning frequency, sampling frequency, and number of data points in
%     the transform (includes negative frequencies).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSNPTS and are restored when the IDFT command is performed.
%
%    Note:
%     - SAC (and thus SAClab for sanity) calculates amp/real/imag spectral 
%       data according to Parseval's theorem.  Dividing records (except for
%       phase!) by npts*delta/2 gives more interesting results for
%       sinusoids.
%
%    Usage: data=dft(data,format,pow2pad)
%
%    Examples:
%
%    See also: idft

% check nargin
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'x'))

% defaults
if(nargin<3 || isempty(pow2pad)); pow2pad=1; end
if(nargin<2 || isempty(format)); format='amph'; end

% check inputs
if(~any(strcmpi(format,{'amph' 'rlim'})))
    error('SAClab:dft:badInput','bad FORMAT string: %s',format)
end
if(~isnumeric(pow2pad) || ~isscalar(pow2pad) || fix(pow2pad)~=pow2pad)
    error('SAClab:dft:badInput','POW2PAD must be a scalar integer')
end

% retreive header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
[b,delta]=gh(data,'b','delta');

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('SAClab:dft:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Time Series File')...
        & ~strcmp(iftype,'General X vs Y file')))
    error('SAClab:dft:illegalOperation',...
        'illegal operation on spectral file')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of datum/components
    [len,ncmp]=size(data(i).x);
    
    % get frequency info
    nspts=2^(nextpow2(len)+pow2pad);
    sb=0; se=1/(delta(i)*2); sdelta=2*se/nspts;
    
    % fft
    data(i).x=delta(i)*fft(data(i).x,nspts);    % SAC compatible
    %data(i).x=2*fft(data(i).x,nspts)/len;      % Accurate sinusoid amplitudes
    
    % expand data to make room for split
    data(i).x(:,(1:ncmp)*2)=data(i).x;
    
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
