function [data]=getspectralcmp(data,cmp)
%GETSPECTRALCMP    Returns the indicated portion of spectral records
%
%    Usage:    data=getspectralcmp(data,'am'|'ph'|'rl'|'im'|'pw'|'cmplx')
%
%    Description:
%     GETSPECTRALCMP(DATA,CMP) extracts the spectral component indicated by
%     CMP.  CMP must be one of the following: 'AM', 'PH', 'RL', 'IM', 'PW'
%     or 'CMPLX'.  Those stand for the amplitude, phase, real, imaginary,
%     power and complex spectras.  CMP may be a list of components as long
%     as there is exactly one entry per record.  Filetype is changed to
%     General X vs Y.  Note that this drops the negative frequencies as
%     they are redundant for real-valued time-series data.  The amplitude
%     spectrum is accordingly doubled.
%
%    Notes:
%     - Using 'CMPLX' will return a complex array.  This will definitely
%       cause issues with other functions.  Use at your own risk.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, IFTYPE, NPTS
%
%    Examples:
%     % This is SAC's KEEPAM:
%     data=getspectralcmp(data,'am');
%
%    See also: KEEPAM, KEEPPH, KEEPRL, KEEPIM, KEEPPW,
%              SPLITRECORDS, DFT,IDFT

%     Version History:
%        June 25, 2009 - initial version
%        Oct. 13, 2009 - does fftshift and adjusts B
%        Oct. 14, 2009 - just return positive freqs, no B or E adjust
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        Aug. 15, 2010 - nargchk fix
%        Feb. 11, 2011 - mass seizmocheck fix
%        Dec. 21, 2011 - add power spectra option, better checkheader usage
%        Feb.  6, 2012 - power spectra output is in dBs
%        Mar. 13, 2012 - use getheader improvements
%        June  1, 2012 - bugfix: multiply amp spectra by 2
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2012 at 15:05 GMT

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
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % valid cmp
    valid.CMP={'am' 'ph' 'rl' 'im' 'pw' 'cmplx'};

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
    [npts,ncmp,sdelta,nspts,iftype]=getheader(data,...
        'npts','ncmp','sdelta','nspts','iftype id');
    npts=npts/2+1; % new npts

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
                    data(i).dep=2*abs(complex(...
                        data(i).dep(:,1:2:end),data(i).dep(:,2:2:end)));
                else
                    data(i).dep=2*data(i).dep(:,1:2:end);
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
            case 'pw'
                % divide by delta factor (parseval's) b/c of squaring
                if(isrlim(i))
                    data(i).dep=complex(data(i).dep(:,1:2:end),...
                        data(i).dep(:,2:2:end));
                else
                    data(i).dep=data(i).dep(:,1:2:end)...
                        .*exp(1j*data(i).dep(:,2:2:end));
                end
                data(i).dep=2*data(i).dep.*conj(data(i).dep)/nspts(i);
                data(i).dep=data(i).dep/sdelta(i);
                
                % convert to dBs
                data(i).dep=10*log10(data(i).dep);
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
        depmen(i)=nanmean(data(i).dep(:));
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
