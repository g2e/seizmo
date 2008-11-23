function [data]=rlim2amph(data)
%RLIM2AMPH    Convert SEIZMO spectral records from RLIM to AMPH
%
%    Description: RLIM2AMPH(DATA) converts SEIZMO real-imaginary records 
%     to amplitude-phase records.  This is particularly useful for 
%     switching between the formats when performing basic operations on 
%     spectral records would require separating the amplitude and phase
%     components.  Records in DATA must be of the spectral variety.  
%     Amplitude-phase records are not altered.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: IFTYPE, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=rlim2amph(data)
%
%    Examples:
%     To simply multiply two records in the frequency domain, they must be
%     converted to real-imaginary first:
%      data=amph2rlim(data)
%      data=multiplyrecords(data(1),data(2))
%      data=rlim2amph(data)
%
%    See also: amph2rlim, dft, idft

%     Version History:
%        June 11, 2008 - initial version
%        July 19, 2008 - removed option, single call to ch, dataless
%                        support, updates DEP* fields, .dep rather than .x,
%                        doc update
%        Oct.  7, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 06:30 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% retreive header info
iftype=getenumdesc(data,'iftype');

% records must be spectral
if(any(~strcmpi(iftype,'Spectral File-Real/Imag')...
        & ~strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:rlim2amph:illegalOperation',...
        'Illegal operation on non-spectral file!');
end

% loop through records
nrecs=numel(data);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % convert or message
    if(strcmpi(iftype(i),'Spectral File-Real/Imag'))
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        temp=complex(data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
        data(i).dep(:,1:2:end)=abs(temp);
        data(i).dep(:,2:2:end)=angle(temp);
        data(i).dep=oclass(data(i).dep);
    end
    
    % dep*
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update filetype
data=changeheader(data,'iftype','Spectral File-Ampl/Phase',...
    'depmax',depmax,'depmin',depmin,'depmen',depmen);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
