function [data]=idft(data)
%IDFT    Performs an inverse discrete fourier transform on SAClab records
%
%    Description: IDFT(DATA) converts SAClab records from the frequency 
%     domain to the time domain using the inverse fast fourier transform.  
%     Output filetype is a Time Series File.
%
%    Header Changes: B, SB, E, DELTA, SDELTA, NPTS, NSPTS
%     In the frequency domain B, DELTA, and NPTS are changed to the 
%     beginning frequency, sampling frequency, and number of data points in
%     the transform (includes negative frequencies).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSNPTS and are restored when this command is performed.
%
%    Notes:
%     - Assumes amplitudes for spectral files follow the convention from 
%       SAC/SAClab's DFT routine.  There is a conversion factor of 
%       npts*delta/2 between this and spectral amplitudes that are more
%       interesting for sinusoids.
%
%    Usage: data=idft(data)
%
%    Examples:
%     To take the derivative of a time-series in the frequency domain:
%      data=idft(mulomega(dft(data)))
%
%    See also: dft, amph2rlim, rlim2amph, divomega, mulomega

%     Version History:
%        ????????????? - Initial Version
%        June 11, 2008 - Added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2008 at 19:50 GMT

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% retreive header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
[b,delta,sb,sdelta,npts,nspts]=...
    gh(data,'b','delta','sb','sdelta','npts','nspts');
e=sb+(nspts-1).*sdelta;

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('SAClab:idft:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Spectral File-Real/Imag')...
        & ~strcmp(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:idft:illegalOperation',...
        'illegal operation on non-spectral file')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % turn back into time domain
    if(strcmp(iftype(i),'Spectral File-Real/Imag'))
        data(i).x=1/sdelta(i)*ifft(complex(data(i).x(:,1:2:end),data(i).x(:,2:2:end)),'symmetric');
        data(i)=ch(data(i),'iftype','Time Series File');
    else
        data(i).x=1/sdelta(i)*ifft(data(i).x(:,1:2:end).*exp(j*data(i).x(:,2:2:end)),'symmetric');
        data(i)=ch(data(i),'iftype','Time Series File');
    end
    
    % truncate to original length and change class back
    data(i).x=oclass(data(i).x(1:nspts(i),:));
end

% update header (note there is no field 'se')
data=ch(data,'b',sb,'e',e,'delta',sdelta,'sb',b,'sdelta',delta,...
    'nspts',npts,'npts',nspts);

end
