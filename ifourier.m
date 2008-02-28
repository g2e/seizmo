function [data]=ifourier(data)
%IFOURIER    Converts spectral seislab records to the time domain
%
%    Description: Converts seislab records from the frequency domain to the
%     time domain using the inverse fast fourier transform.  Output file
%     type is timeseries file.
%
%    Note:
%     - Assumes amplitudes for amplitude-phase spectral files follow the
%       convention from seislab's fourier routine.  This adds a factor of
%       1/npts to frequency amplitudes compared to SAC.
%
%    Usage: data=ifourier(data)
%
%    See also: fourier

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% retreive header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
[b,delta,sb,sdelta,npts,nspts]=...
    gh(data,'b','delta','sb','sdelta','npts','nspts');
e=sb+(nspts-1)*sdelta;

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('seislab:fourier:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Spectral File-Real/Imag'))...
        && any(~strcmp(iftype,'Spectral File-Ampl/Phase')))
    error('seislab:fourier:illegalOperation',...
        'illegal filetype for this operation')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % turn back into time domain
    if(strcmp(iftype(i),'Spectral File-Real/Imag'))
        data(i).x=nspts(i)*real(ifft(complex(data(i).x(:,1),data(i).x(:,2))));
        data(i)=ch(data(i),'iftype','itime');
    else
        data(i).x=nspts(i)/2*real(ifft(data(i).x(:,1).*exp(j*data(i).x(:,2))));
        data(i)=ch(data(i),'iftype','itime');
    end
    
    % truncate to original length and change class back
    data(i).x=oclass(data(i).x(1:nspts(i),:));
end

% update header (note there is no field 'se')
data=ch(data,'b',sb,'e',e,'delta',sdelta,'sb',b,'sdelta',delta,...
    'nspts',npts,'npts',nspts);

end
