function [data]=divomega(data)
%DIVOMEGA    Perform integration in the frequency domain on SAClab records
%
%    Description: DIVOMEGA(DATA) divides each point of a spectral file by
%     its frequency given by:
%       OMEGA=2.0 * PI * FREQ
%     This is analogous to integrating the equivalent time series file.
%     The spectral file can be in either amplitude-phase or real-imaginary
%     format.
%     
%     This is often convenient with normal data but is critical when
%     obtaining the DFT of data whose spectra ranges over many orders of
%     magnitude.  For example suppose you have prewhitened a data file by
%     using the DIF command, and then taken its transform using the DFT
%     command.  The effect of the differentiation in the time domain can be
%     removed by an integration in the frequency domain using this command.
%
%    Header Changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  data=divomega(data)
%
%    Examples:
%     Integrate time series data in the frequency domain vs time domain:
%      data=idft(divomega(dft(data)))
%      data=integrt(data)
%
%    See also: mulomega, dft, idft

%     Version History:
%        ????????????? - Initial Version
%        June 11, 2008 - Documentation Cleanup
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
[e,delta,npts]=gh(data,'e','delta','npts');
npts2=npts/2;
npts21=npts2-1;

% check leven,iftype
if(any(~strcmp(leven,'true')))
    error('SAClab:divomega:illegalOperation',...
        'illegal operation on unevenly spaced record')
elseif(any(~strcmp(iftype,'Spectral File-Real/Imag')...
        & ~strcmp(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:divomega:illegalOperation',...
        'illegal operation on a non-spectral file')
end

% loop through records
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % integrate
    cols=size(data(i).x,2);
    omega=[0 1./(2*pi*[linspace(delta(i),e(i),npts2(i)) ...
        linspace(e(i)-delta(i),delta(i),npts21(i))])].';
    if(strcmp(iftype(i),'Spectral File-Real/Imag'))
        % rlim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0Hz real/imag == 0 else
        % real=imag/omega & imag=-real/omega
        data(i).x(:,[1:2:end 2:2:end])=oclass(...
            [data(i).x(:,2:2:end).*omega(:,ones(1,cols/2))...
            -data(i).x(:,1:2:end).*omega(:,ones(1,cols/2))]);
    else
        % amph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % amp at 0Hz == 0 else amp=amp/omega
        % phase=phase+pi/2 at -Hz and 0Hz else phase=phase-pi/2
        data(i).x(:,1:2:end)=data(i).x(:,1:2:end).*omega(:,ones(1,cols/2));
        data(i).x(2:(npts2+1),2:2:end)=data(i).x(2:(npts2+1),2:2:end)-pi/2;
        data(i).x([1 npts2+2:end],2:2:end)=data(i).x([1 npts2+2:end],2:2:end)+pi/2;
        data(i).x=oclass(data(i).x);
    end
    
    % update header
    data(i)=ch(data(i),'depmax',max(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmen',mean(data(i).x(:)));
end

end
