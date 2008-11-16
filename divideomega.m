function [data]=divomega(data)
%DIVOMEGA    Perform integration in the frequency domain on SEIZMO records
%
%    Description: DIVOMEGA(DATA) basically divides each point in the 
%     dependent component(s) of spectral files by:
%       OMEGA=2.0 * PI * FREQ
%     to perform the equivalent of integration in the time domain.  This is
%     particularly handy when working with spectral data as it avoids
%     the forward and inverse fourier transform necessary for time domain 
%     integration.  It is also useful for reducing the dynamic range of 
%     spectral data.
%
%    Notes:
%     - Read the source code below for a better description of the
%       operations performed for frequency-domain integration.
%
%    System requirements: Matlab 7
%
%    Header Changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=divomega(data)
%
%    Examples:
%     Integrate spectral data in the time domain vs frequency domain:
%      data=dft(integrt(idft(data)))
%      data=divomega(data)
%
%    See also: mulomega, dft, idft

%     Version History:
%        May  12, 2008 - initial version
%        June 11, 2008 - doc cleanup
%        July  8, 2008 - doc update, single ch call, .dep rather than .x
%        July 19, 2008 - doc update, dataless support
%        Oct.  7, 2008 - minor code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 02:10 GMT

% todo:

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
if(any(~strcmpi(leven,'true')))
    error('seizmo:divomega:illegalOperation',...
        'Illegal operation on unevenly spaced record!');
elseif(any(~strcmpi(iftype,'Spectral File-Real/Imag')...
        & ~strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:divomega:illegalOperation',...
        'Illegal operation on a non-spectral file!');
end

% number of records
nrecs=numel(data);

% loop through records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
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
