function [data]=getspectralcmp(data,cmp)
%GETSPECTRALCMP    Returns the indicated portion of spectral records
%
%    Usage:    data=getspectralcmp(data,'am'|'ph'|'rl'|'im'|'cmplx')
%
%    Description: GETSPECTRALCMP(DATA,CMP) extracts the spectral component
%     indicated by CMP.  CMP must be one of the following: 'AM', 'PH',
%     'RL', 'IM', or 'CMPLX'.  CMP may be a list of components as long as
%     there is exactly one entry per record.  Filetype is changed to
%     General X vs Y.
%
%    Notes:
%     - Using 'CMPLX' will return a complex array.  This will definitely
%       cause issues with other functions.  Use at your own risk.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE
%
%    Examples:
%     This is SAC's KEEPAM:
%      data=getspectralcmp(data,'am');
%
%    See also: KEEPAM, KEEPPH, KEEPRL, KEEPIM, SPLITRECORDS, DFT, IDFT

%     Version History:
%        June 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:05 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% number of records
nrecs=numel(data);

% valid cmp
valid.CMP={'am' 'ph' 'rl' 'im' 'cmplx'};

% check/prepare cmp
if(iscellstr(cmp))
    cmp=char(cmp);
end
if(isempty(cmp) || ~ischar(cmp) || ~any(size(cmp,1)==[1 nrecs]) ...
        || ~isempty(setdiff(lower(cmp),valid.CMP)))
    error('seizmo:dft:badInput',...
        ['CMP must be one of the following:\n' ...
        sprintf('%s ',valid.CMP{:})]);
end
if(size(cmp,1)==1)
    cmp=cmp(ones(nrecs,1),:);
end
cmp=cellstr(lower(cmp));

% get filetype
iftype=getenumdesc(data,'iftype');
ncmp=getncmp(data);

% require spectral records
if(any(~strcmpi(iftype,'Spectral File-Real/Imag')...
        & ~strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:idft:illegalOperation',...
        'Illegal operation on non-spectral file!');
end

% logical array for filetype
isrlim=strcmpi(iftype,'Spectral File-Real/Imag');

% loop over records
iftype='ixy';
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless (but decrease columns)
    if(isempty(data(i).dep)); data(i).dep=zeros(0,ncmp); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % which component
    switch cmp{i}
        case 'am'
            if(isrlim(i))
                data(i).dep=abs(complex(...
                    data(i).dep(:,1:2:end),data(i).dep(:,2:2:end)));
            else
                data(i).dep=data(i).dep(:,1:2:end);
            end
        case 'ph'
            if(isrlim(i))
                data(i).dep=angle(complex(...
                    data(i).dep(:,1:2:end),data(i).dep(:,2:2:end)));
            else
                data(i).dep=data(i).dep(:,2:2:end);
            end
        case 'rl'
            if(isrlim(i))
                data(i).dep=data(i).dep(:,1:2:end);
            else
                data(i).dep=real(...
                    data(i).dep(:,1:2:end).*exp(j*data(i).dep(:,2:2:end)));
            end
        case 'im'
            if(isrlim(i))
                data(i).dep=data(i).dep(:,2:2:end);
            else
                data(i).dep=imag(...
                    data(i).dep(:,1:2:end).*exp(j*data(i).dep(:,2:2:end)));
            end
        case 'cmplx'
            if(isrlim(i))
                data(i).dep=complex(...
                    data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
            else
                data(i).dep=...
                    data(i).dep(:,1:2:end).*exp(j*data(i).dep(:,2:2:end));
            end
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'iftype',iftype);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end

