function [data]=amph2rlim(data,option,ignore)
%AMPH2RLIM    Convert SAClab spectral records from AMPH to RLIM
%
%    Description: AMPH2RLIM(DATA) converts SAClab amplitude-phase records 
%     to real-imaginary records.  This is particularly useful when
%     performing basic operations on spectral records which would otherwise
%     require separating the amplitude and phase components.  Records in 
%     DATA must be of the spectral variety.  Real-imaginary records are
%     not altered (a warning message is given unless option 
%     'ignore_preconverted' is set to TRUE).
%
%    Notes:
%
%    Requirements: Matlab 7
%
%    Supported types: ALL
%
%    Header changes: IFTYPE
%
%    Usage:  data=amph2rlim(data)
%            data=amph2rlim(data,'ignore_preconverted',true|false)
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
%        June 11, 2008 - Initial Version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2008 at 20:30 GMT

% check nargin
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'x'))

% check option
if(nargin==1)
    ignore=false;
elseif(nargin==2 || ~strcmp(option,'ignore_preconverted'))
    error('SAClab:amph2rlim:badOption','Bad option: %s',option)
end 

% retreive header info
iftype=genumdesc(data,'iftype');

% records must be spectral
if(any(~strcmp(iftype,'Spectral File-Real/Imag')...
        & ~strcmp(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:amph2rlim:illegalOperation',...
        'illegal operation on non-spectral file')
end

% loop through records
for i=1:length(data)
    % convert or message
    if(strcmp(iftype(i),'Spectral File-Ampl/Phase'))
        oclass=str2func(class(data(i).x));
        data(i).x=double(data(i).x);
        temp=data(i).x(:,1:2:end).*exp(j*data(i).x(:,2:2:end));
        data(i).x(:,1:2:end)=real(temp);
        data(i).x(:,2:2:end)=imag(temp);
        data(i).x=oclass(data(i).x);
    else
        if(~ignore)
            warning('SAClab:amph2rlim:unnecessaryOperation',...
            'Unnecessary operation on real-imaginary file %s',data(i).name)
        end
    end
    
    % update filetype
    data=ch(data,'iftype','Spectral File-Real/Imag');
end

end
