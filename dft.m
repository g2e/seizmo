function [data]=dft(data,format,pow2pad)
%DFT    Performs a discrete fourier transform on SAClab data records
%
%    Description: DFT(DATA,FORMAT) converts SAClab records from the time 
%     domain to the frequency domain using a discrete fourier transform
%     (Matlab's fft).  Following SAC formatting, an option FORMAT can be
%     given to select between real-imaginary and amplitude-phase formats
%     (FORMAT can be either 'amph' or 'rlim' - the default is 'amph').  
%
%     DFT(DATA,FORMAT,POW2PAD) lets the power of 2 zero-padding be adjusted
%     using an integer POW2PAD according to the formula:
%                fftlength=2^(nextpow2(NPTS)+POW2PAD)
%     The default value is 1.
%
%    Notes:
%     - SAC (and thus SAClab for sanity) calculates spectral data according
%       to Parseval's theorem.  Dividing records (except for phase!) by 
%       npts*delta/2 gives results that may better match the amplitudes
%       of sinusoids.
%
%    System requirements: Matlab 7
%
%    Data requirements: Evenly Spaced; Time Series or General X vs Y
%
%    Header Changes: B, SB, E, DELTA, SDELTA, NPTS, NSPTS
%                    DEPMEN, DEPMIN, DEPMAX
%     In the frequency domain B, DELTA, and NPTS are changed to the 
%     beginning frequency, sampling frequency, and number of data points in
%     the transform (includes negative frequencies).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSNPTS and are restored when the IDFT command is performed.
%
%    Usage: data=dft(data)
%           data=dft(data,format)
%           data=dft(data,format,pow2pad)
%
%    Examples:
%     To take the derivative of a time-series in the frequency domain:
%      data=idft(mulomega(dft(data)))
%
%    See also: idft, amph2rlim, rlim2amph, divomega, mulomega

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 28, 2008 - seischk support and code cleanup
%        Mar.  2, 2008 - readibility change
%        Mar.  3, 2008 - now compatible with SAC
%        Mar.  4, 2008 - cleaned up errors and warnings
%        May  12, 2008 - name changed to dft
%        June 11, 2008 - added example
%        June 29, 2008 - documentation update, .dep rather than .x
%        July 18, 2008 - dataless support, one ch call, updates DEP* fields
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 18, 2008 at 23:15 GMT

% todo:
% 

% check nargin
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'dep'))

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
if(any(~strcmpi(leven,'true')))
    error('SAClab:dft:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('SAClab:dft:illegalOperation',...
        'illegal operation on spectral/xyz record')
end

% output type
if(strcmpi(format,'rlim')); iftype='Spectral File-Real/Imag';
else iftype='Spectral File-Ampl/Phase';
end

% loop through records
sb=nan(nrecs,1); se=sb; sdelta=sb; nspts=sb; len=sb;
depmen=sb; depmin=sb; depmax=sb;
for i=1:numel(data)
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % number of datum/components
    [len(i),ncmp]=size(data(i).dep);
    
    % get frequency info
    nspts=2^(nextpow2(len(i))+pow2pad);
    sb(i)=0; se(i)=1/(delta(i)*2); sdelta=2*se(i)/nspts(i);
    
    % fft
    data(i).dep=delta(i)*fft(data(i).dep,nspts(i));    % SAC compatible
    %data(i).dep=2*fft(data(i).dep,nspts(i))/len(i);   % better sinusoid amplitudes
    
    % expand data to make room for split
    data(i).dep(:,(1:ncmp)*2)=data(i).dep;
    
    % split complex by desired filetype
    if(strcmpi(format,'rlim'))
        data(i).dep(:,1:2:end)=real(data(i).dep(:,2:2:end));
        data(i).dep(:,2:2:end)=imag(data(i).dep(:,2:2:end));
    else
        data(i).dep(:,1:2:end)=abs(data(i).dep(:,2:2:end));
        data(i).dep(:,2:2:end)=angle(data(i).dep(:,2:2:end));
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header (note there is no field 'se')
data=ch(data,'b',sb,'e',se,'delta',sdelta,'sb',b,...
    'sdelta',delta,'nspts',len,'npts',nspts,'iftype',iftype,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
