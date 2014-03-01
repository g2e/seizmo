function [data]=removepolynomial(data,degree)
%REMOVEPOLYNOMIAL    Remove polynomial fit from SEIZMO records
%
%    Usage:    data=removepolynomial(data,degree)
%
%    Description:
%     DATA=REMOVEPOLYNOMIAL(DATA,DEGREE) removes the polynomial fit with
%     polynomial degree DEGREE from SEIZMO records in DATA.  For multi-
%     component records, each component is dealt with separately.  DEGREE
%     may be a scalar or an array of size NRECSxNCMP to specify different
%     values for each record & component.
%
%    Notes:
%     - Degree 0 == REMOVEMEAN
%       Degree 1 == REMOVETREND
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Check out the difference various degree polynomials make:
%     plot1(removepolynomial(data(ones(1,10)),0:9))
%
%    See also: REMOVESPLINE, REMOVEMEAN, REMOVETREND, GETPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS, GETSPLINE, POLYFIT, POLYVAL

%     Version History:
%        June 24, 2009 - initial version
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%        Jan. 30, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  2, 2010 - fix verbose bug, versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Mar. 13, 2012 - doc update, seizmocheck fix, leven fix, use
%                        getheader improvements
%        Jan. 21, 2014 - minor doc fix
%        Jan. 22, 2014 - bugfix: poly. order was actually degree, changed
%                        all instances of order to degree, allow per cmp
%                        degree specification, turn off polyfit warnings
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt polynomial removal
try
    % turn off the badly conditioned polynomial warnings
    warning('off','MATLAB:polyfit:RepeatedPoints');
    warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
    
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check degree
    if(~isnumeric(degree) || any(degree(:)~=fix(degree(:))) ...
            || (~isscalar(degree) && mod(numel(degree),nrecs)~=0) ...
            || ~any(size(degree,1)==[1 nrecs]))
        error('seizmo:removepolynomial:badOrder',...
            'DEGREE must be a scalar or an array of integers.');
    end
    if(size(degree,2)==nrecs && size(degree,1)==1)
        % make row vector a column vector
        degree=degree(:);
    end
    if(~any(size(degree,1)==[1 nrecs]))
        error('seizmo:removepolynomial:badOrder',...
            'DEGREE must be a scalar or a NRECSxNCMP array.');
    end

    % header info
    [delta,npts,ncmp,leven]=getheader(data,...
        'delta','npts','ncmp','leven lgc');
    leven=~strcmpi(leven,'false');
    
    % expand degree
    if(isscalar(degree))
        degree=degree(ones(nrecs,max(ncmp)));
    elseif(size(degree,2)==1 && max(ncmp)>1)
        degree=degree(:,ones(1,max(ncmp)));
    elseif(size(degree,2)==max(ncmp) && size(degree,1)==1)
        degree=degree(ones(nrecs,1),:);
    end
    
    % detail message
    if(verbose)
        disp('Removing Polynomial Fit from Record(s)');
        print_time_left(0,nrecs);
    end

    % remove polynomial and update header
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
        if(leven(i))
            time=(0:delta(i):delta(i)*(npts(i)-1)).';
            for j=1:ncmp(i)
                data(i).dep(:,j)=data(i).dep(:,j) ...
                    -polyval(polyfit(time,data(i).dep(:,j),degree(i,j)),...
                    time);
            end
        % unevenly spaced
        else
            for j=1:ncmp(i)
                data(i).dep(:,j)=data(i).dep(:,j)...
                    -polyval(polyfit(double(data(i).ind),...
                    data(i).dep(:,j),degree(i,j)),double(data(i).ind));
            end
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % adjust header
        depmen(i)=nanmean(data(i).dep(:));
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
    
    % turn on the badly conditioned polynomial warnings
    warning('on','MATLAB:polyfit:RepeatedPoints');
    warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % turn on the badly conditioned polynomial warnings
    warning('on','MATLAB:polyfit:RepeatedPoints');
    warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
    
    % rethrow error
    error(lasterror);
end

end
