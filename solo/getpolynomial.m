function [p]=getpolynomial(data,order)
%GETPOLYNOMIAL    Get polynomial fit to SEIZMO records
%
%    Usage:    p=getpolynomial(data,order)
%
%    Description:
%     P=GETPOLYNOMIAL(DATA,ORDER) gets the polynomial trend of order ORDER
%     from SEIZMO records.  For multi-component records, each component is
%     dealt with separately.  Polynomial coefficients are returned in cell
%     array P, with each cell corresponding to each record in DATA and each
%     row in each cell giving the coefficients for each component.
%
%    Notes:
%
%    Examples:
%     % Get various polynomial fits to a record:
%     getpolynomial(data(ones(1,5)),0:4)
%
%    See also: REMOVEMEAN, REMOVETREND, REMOVEPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS

%     Version History:
%        June 24, 2009 - initial version
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        Feb.  2, 2010 - versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Mar. 13, 2012 - doc update, seizmocheck fix, leven fix, use
%                        getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check input
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt polynomial fit
try
    % check headers
    data=checkheader(data);

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check order
    if(~isnumeric(order) || any(order~=fix(order)) ...
            || ~any(numel(order)==[1 nrecs]))
        error('seizmo:getpolynomial:badOrder',...
            'ORDER must be a scalar or an array of integers.');
    end
    if(isscalar(order))
        order(1:nrecs,1)=order;
    end

    % header info
    [delta,npts,leven]=getheader(data,'delta','npts','leven lgc');
    leven=~strcmpi(leven,'false');
    
    % detail message
    if(verbose)
        disp('Getting Polynomial Fit to Record(s)');
        print_time_left(0,nrecs);
    end

    % remove trend and update header
    p=cell(nrecs,1);
    for i=1:numel(data)
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % convert to double precision
        data(i).dep=double(data(i).dep);

        % evenly spaced
        if(leven(i))
            time=((0:npts(i)-1)*delta(i)).';
            for j=1:size(data(i).dep,2)
                p{i}(j,1:order(i)+1)=...
                    polyfit(time,data(i).dep(:,j),order(i));
            end
        % unevenly spaced
        else
            for j=1:size(data(i).dep,2)
                p{i}(j,1:order(i)+1)=...
                    polyfit(double(data(i).ind),data(i).dep(:,j),order(i));
            end
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
