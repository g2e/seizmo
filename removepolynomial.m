function [data]=removepolynomial(data,order)
%REMOVEPOLYNOMIAL    Remove polynomial trend from SEIZMO records
%
%    Usage:    data=removepolynomial(data,order)
%
%    Description: REMOVEPOLYNOMIAL(DATA,ORDER) removes the polynomial trend
%     of order ORDER from SEIZMO records.  For multi-component
%     records, each component is dealt with separately.  It is highly
%     recommended to combine this command with filtering operations to
%     reduce edge effects that may lead to poor data quality.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Check out the difference various order polynomials make:
%      plot1(removepolynomial(data(ones(1,10)),1:10))
%
%    See also: REMOVEMEAN, REMOVETREND, GETPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS

%     Version History:
%        June 24, 2009 - initial version
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%        Jan. 30, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  2, 2010 - fix verbose bug, versioninfo caching
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  2, 2010 at 21:25 GMT

% todo:

% check input
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt polynomial removal
try
    % check headers (versioninfo cache update)
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check order
    if(~isnumeric(order) || any(order~=fix(order)) ...
            || ~any(numel(order)==[1 nrecs]))
        error('seizmo:removepolynomial:badOrder',...
            'ORDER must be a scalar or an array of integers.');
    end
    if(isscalar(order))
        order(1:nrecs,1)=order;
    end

    % header info
    [delta,npts,ncmp]=getheader(data,'delta','npts','ncmp');
    leven=getlgc(data,'leven');
    
    % detail message
    if(verbose)
        disp('Removing Polynomial from Record(s)');
        print_time_left(0,nrecs);
    end

    % remove trend and update header
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:numel(data)
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % evenly spaced
        if(strcmp(leven(i),'true'))
            time=((0:npts(i)-1)*delta(i)).';
            for j=1:ncmp(i)
                data(i).dep(:,j)=data(i).dep(:,j) ...
                    -polyval(polyfit(time,data(i).dep(:,j),order(i)),time);
            end
            % unevenly spaced
        else
            for j=1:ncmp(i)
                data(i).dep(:,j)=data(i).dep(:,j)...
                    -polyval(polyfit(double(data(i).ind),...
                    data(i).dep(:,j),order(i)),double(data(i).ind));
            end
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % adjust header
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % adjust header
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
