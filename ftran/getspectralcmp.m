function [data]=getspectralcmp(data,cmp)
%GETSPECTRALCMP    Returns the indicated portion of spectral records
%
%    Usage:    data=getspectralcmp(data,'am'|'ph'|'rl'|'im'|'cmplx')
%
%    Description: GETSPECTRALCMP(DATA,CMP) extracts the spectral component
%     indicated by CMP.  CMP must be one of the following: 'AM', 'PH',
%     'RL', 'IM', or 'CMPLX'.  CMP may be a list of components as long as
%     there is exactly one entry per record.  Filetype is changed to
%     General X vs Y.  Note that this drops the negative frequencies.
%
%    Notes:
%     - Using 'CMPLX' will return a complex array.  This will definitely
%       cause issues with other functions.  Use at your own risk.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE, NPTS
%
%    Examples:
%     This is SAC's KEEPAM:
%      data=getspectralcmp(data,'am');
%
%    See also: KEEPAM, KEEPPH, KEEPRL, KEEPIM, SPLITRECORDS, DFT, IDFT

%     Version History:
%        June 25, 2009 - initial version
%        Oct. 13, 2009 - does fftshift and adjusts B
%        Oct. 14, 2009 - just return positive freqs, no B or E adjust
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        Aug. 15, 2010 - nargchk fix
%        Feb. 11, 2011 - mass seizmocheck fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt cmp extraction
try
    % check headers
    data=checkheader(data);

    % verbosity
    verbose=seizmoverbose;

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
        error('seizmo:getspectralcmp:badInput',...
            ['CMP must be one of the following:\n' ...
            sprintf('%s ',valid.CMP{:})]);
    end
    if(size(cmp,1)==1)
        cmp=cmp(ones(nrecs,1),:);
    end
    cmp=cellstr(lower(cmp));

    % get header info
    iftype=getenumid(data,'iftype');
    [npts,ncmp]=getheader(data,'npts','ncmp');
    npts=npts/2+1; % new npts

    % require spectral records
    if(any(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph')))
        error('seizmo:getspectralcmp:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'irlim') & ~strcmpi(iftype,'iamph'))) ...
            '\nDatatype of record(s) in DATA must be spectral!']);
    end

    % logical array for filetype
    isrlim=strcmpi(iftype,'irlim');
    
    % detail message
    if(verbose)
        disp('Extracting Component from Spectral Record(s)');
        print_time_left(0,nrecs);
    end

    % loop over records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless (but decrease columns)
        if(isempty(data(i).dep))
            data(i).dep=zeros(0,ncmp(i));
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

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
                    data(i).dep=real(data(i).dep(:,1:2:end)...
                        .*exp(1j*data(i).dep(:,2:2:end)));
                end
            case 'im'
                if(isrlim(i))
                    data(i).dep=data(i).dep(:,2:2:end);
                else
                    data(i).dep=imag(data(i).dep(:,1:2:end)...
                        .*exp(1j*data(i).dep(:,2:2:end)));
                end
            case 'cmplx'
                if(isrlim(i))
                    data(i).dep=complex(...
                        data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
                else
                    data(i).dep=data(i).dep(:,1:2:end)...
                        .*exp(1j*data(i).dep(:,2:2:end));
                end
        end

        % change class back
        data(i).dep=oclass(data(i).dep(1:npts(i),:));

        % dep*
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,'iftype','ixy','npts',npts,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
