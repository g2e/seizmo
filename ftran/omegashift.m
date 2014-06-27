function [data]=omegashift(data,shift)
%OMEGASHIFT    Applies time shift to SEIZMO data in the frequency domain
%
%    Usage:    data=omegashift(data,shift)
%
%    Description:
%     DATA=OMEGASHIFT(DATA,SHIFT) shifts the data in SEIZMO struct DATA by
%     SHIFT seconds without changing the timing.  A positive SHIFT will
%     delay the data while a negative SHIFT will advance the data.  Since
%     the operation is done in the frequency domain there is a circular
%     shift occurring (that is data from the end wraps to the front and
%     data from the front wraps to the end) but it is complicated by
%     records being zero-padded.  SHIFT may be a scalar or a vector of
%     values (1 per record in DATA).
%
%    Notes:
%     - TIMESHIFT & OMEGASHIFT are completely different operations!  The
%       timing of the data is not preserved in OMEGASHIFT (you may want to
%       use TIMESHIFT to do so).
%
%    Examples:
%     % Stack P wave arrivals:
%     b=getheader(data,'b');
%     P=findpicks(arrivals2picks(data,'P,Pdiff'),'P,Pdiff',1);
%     plot1(addrecords(idft(...
%         omegashift(dft(data,'rlim'),b-P)),'npts','pad'))
%
%    See also: INTERPOLATE, TIMESHIFT, OMEGADIVIDE, OMEGAMULTIPLY, DFT,
%              IDFT, OMEGAANALYTIC, OMEGAGAUSSIAN

%     Version History:
%        Feb.  4, 2012 - initial version
%        Mar. 15, 2012 - fix example
%        May  31, 2012 - minor doc update
%        Feb. 14, 2013 - use strcmpi for consistency
%        Mar.  1, 2014 - maketime fix
%        June 26, 2014 - bugfix: shift==0 no longer causes error, irlim
%                        shift no longer errors due to line wrapping issues
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt frequency domain time shift
try
    % check headers
    data=checkheader(data,...
        'NONSPECTRAL_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check shift
    if(~isnumeric(shift) || ~isreal(shift))
        error('seizmo:omegashift:badInput',...
            'SHIFT must be a real-value in seconds!');
    elseif(~any(numel(shift)~=[1 nrecs]))
        error('seizmo:omegashift:badInput',...
            'SHIFT must be equal in size to DATA or scalar!');
    end
    if(isscalar(shift)); shift(1:nrecs)=shift; end
    
    % retrieve header info
    [delta,npts,iftype]=getheader(data,'delta','npts','iftype id');
    
    % detail message
    if(verbose)
        disp('Time Shifting Record(s) in the Frequency Domain');
        print_time_left(0,nrecs);
    end
    
    % loop through records
    [depmin,depmen,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % time shift
        cols=size(data(i).dep,2)/2;
        wt=2*pi*delta(i);
        wt=shift(i)*(0:wt:wt*(npts(i)-1)).';
        if(strcmpi(iftype(i),'irlim'))
            data(i).dep(:,[1:2:end 2:2:end])=...
                [data(i).dep(:,1:2:end).*cos(wt(:,ones(1,cols)))...
                + data(i).dep(:,2:2:end).*sin(wt(:,ones(1,cols))) ...
                data(i).dep(:,2:2:end).*cos(wt(:,ones(1,cols)))...
                - data(i).dep(:,1:2:end).*sin(wt(:,ones(1,cols)))];
        else % iamph
            data(i).dep(:,2:2:end)=...
                data(i).dep(:,2:2:end)-mod(wt(:,ones(1,cols)),2*pi);
        end

        % change class back
        data(i).dep=oclass(data(i).dep);
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
