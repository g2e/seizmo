function [data]=rlim2amph(data)
%RLIM2AMPH    Convert SAClab spectral records from RLIM to AMPH
%
%    Description: RLIM2AMPH(DATA) converts SAClab real-imaginary records 
%     to amplitude-phase records.  This is particularly useful for 
%     switching between the formats when performing basic operations on 
%     spectral records which would otherwise require separating the
%     amplitude and phase components.  Records in DATA must be of the 
%     spectral variety.  Amplitude-phase records are not altered (a warning
%     message is given).
%
%    Header changes: IFTYPE
%
%    Usage:  data=rlim2amph(data)
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
%    See also: amph2rlim, dft, idft

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
    error('SAClab:rlim2amph:illegalOperation',...
        'illegal operation on non-spectral file')
end

% loop through records
for i=1:length(data)
    % convert or message
    if(strcmp(iftype(i),'Spectral File-Real/Imag'))
        oclass=str2func(class(data(i).x));
        data(i).x=double(data(i).x);
        temp=complex(data(i).x(:,1:2:end),data(i).x(:,2:2:end));
        data(i).x(:,1:2:end)=abs(temp);
        data(i).x(:,2:2:end)=angle(temp);
        data(i).x=oclass(data(i).x);
    else
        if(~ignore)
            warning('SAClab:rlim2amph:unnecessaryOperation',...
            'Unnecessary operation on ampltude-phase file %s',data(i).name)
        end
    end
    
    % update filetype
    data=ch(data,'iftype','Spectral File-Ampl/Phase');
end

end
