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
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: Evenly sampled, spectral records
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
%        May  12, 2008 - initial version
%        June 11, 2008 - documentation cleanup
%        July  8, 2008 - documentation update, single ch call
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July  8, 2008 at 17:50 GMT

% todo:
%

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

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

% number of records
nrecs=numel(data);

% loop through records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % integrate
    cols=size(data(i).dep,2);
    omega=[0 1./(2*pi*[linspace(delta(i),e(i),npts2(i)) ...
        linspace(e(i)-delta(i),delta(i),npts21(i))])].';
    if(strcmp(iftype(i),'Spectral File-Real/Imag'))
        % rlim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0Hz real/imag == 0 else
        % real=imag/omega & imag=-real/omega
        data(i).dep(:,[1:2:end 2:2:end])=oclass(...
            [data(i).dep(:,2:2:end).*omega(:,ones(1,cols/2))...
            -data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2))]);
    else
        % amph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % amp at 0Hz == 0 else amp=amp/omega
        % phase=phase+pi/2 at -Hz and 0Hz else phase=phase-pi/2
        data(i).dep(:,1:2:end)=...
            data(i).dep(:,1:2:end).*omega(:,ones(1,cols/2));
        data(i).dep(2:(npts2+1),2:2:end)=...
            data(i).dep(2:(npts2+1),2:2:end)-pi/2;
        data(i).dep([1 npts2+2:end],2:2:end)=...
            data(i).dep([1 npts2+2:end],2:2:end)+pi/2;
        data(i).dep=oclass(data(i).dep);
    end
    
    % dep*
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmax',depmax,'depmin',depmin,'depmen',depmen);

end
