function [data]=amph2rlim(data)
%AMPH2RLIM    Convert SAClab spectral records from AMPH to RLIM
%
%    Description: AMPH2RLIM(DATA) converts SAClab amplitude-phase records 
%     to real-imaginary records.  This is particularly useful when
%     performing basic operations on spectral records which would otherwise
%     require separating the amplitude and phase components.  Records in 
%     DATA must be of the spectral variety.  Real-imaginary records are
%     not altered.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Data requirements: Spectral records only.
%
%    Header changes: IFTYPE, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  data=amph2rlim(data)
%
%    Examples:
%     To simply multiply two records in the frequency domain, they must be
%     converted to real-imaginary first (the operation can be done on 
%     amplitude-phase records but requires working on the components 
%     independently):
%      data=amph2rlim(data)
%      data=mulf(data(1),data(2))
%      data=rlim2amph(data)
%
%    See also: rlim2amph, dft, idft

%     Version History:
%        June 11, 2008 - initial version
%        June 20, 2008 - minor documentation update
%        June 28, 2008 - fixed call to ch, removed option,
%                        documentation update, .dep rather than .x
%        July 19, 2008 - dataless support, updates DEP* fields
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 19, 2008 at 05:30 GMT

% todo:
% 

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% retreive header info
iftype=genumdesc(data,'iftype');

% records must be spectral
if(any(~strcmpi(iftype,'Spectral File-Real/Imag')...
        & ~strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:amph2rlim:illegalOperation',...
        'illegal operation on non-spectral record(s)')
end

% loop through records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:numel(data)
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % convert or message
    if(strcmpi(iftype(i),'Spectral File-Ampl/Phase'))
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        temp=data(i).dep(:,1:2:end).*exp(j*data(i).dep(:,2:2:end));
        data(i).dep(:,1:2:end)=real(temp);
        data(i).dep(:,2:2:end)=imag(temp);
        data(i).dep=oclass(data(i).dep);
    end
    
    % dep*
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update filetype
data=ch(data,'iftype','Spectral File-Real/Imag',...
    'depmax',depmax,'depmin',depmin,'depmen',depmen);

end
