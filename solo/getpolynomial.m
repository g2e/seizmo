function [data]=getpolynomial(data,degree,pflag)
%GETPOLYNOMIAL    Get polynomial fit to SEIZMO records
%
%    Usage:    data=getpolynomial(data,degree)
%              p=getpolynomial(data,degree,'p')
%
%    Description:
%     DATA=GETPOLYNOMIAL(DATA,DEGREE) returns the polynomial trend of
%     degree DEGREE of SEIZMO records in DATA.  For multi-component
%     records, each component is dealt with separately.  DEGREE may be a
%     scalar or an array of size NRECSxNCMP to specify different values for
%     each record & component.
%
%     P=GETPOLYNOMIAL(DATA,DEGREE,'P') exports the polynomial trend of
%     degree DEGREE from SEIZMO records as polynomial coefficients.  For
%     multi-component records, each component is dealt with separately.
%     Polynomial coefficients are returned in a cell array P of size
%     NRECSxNCMP, with each cell corresponding to each record/component in
%     DATA (e.g., P{i,j} gives the coefficients for the jth component of
%     the ith record).
%
%    Notes:
%     - Degree 0 == MEAN
%       Degree 1 == LINEAR TREND
%
%    Examples:
%     % Get various polynomial fits to a record:
%     plot0(getpolynomial(data(ones(1,5)),0:4))
%
%    See also: GETSPLINE, REMOVEMEAN, REMOVETREND, REMOVEPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS, REMOVESPLINE, POLYVAL, POLYFIT

%     Version History:
%        June 24, 2009 - initial version
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging
%        Feb.  2, 2010 - versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Mar. 13, 2012 - doc update, seizmocheck fix, leven fix, use
%                        getheader improvements
%        Jan. 22, 2014 - data output by default, poly. coeff. option, doc
%                        update, turn off polyfit warnings
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% default pflag is empty
if(nargin<3); pflag=[]; end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt polynomial fit
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
        error('seizmo:getpolynomial:badOrder',...
            'DEGREE must be a scalar or an array of integers.');
    end
    if(size(degree,2)==nrecs && size(degree,1)==1)
        % make row vector a column vector
        degree=degree(:);
    end
    if(~any(size(degree,1)==[1 nrecs]))
        error('seizmo:getpolynomial:badOrder',...
            'DEGREE must be a scalar or a NRECSxNCMP array.');
    end
    
    % check pflag
    if(~isempty(pflag) && ischar(pflag) && isequal(lower(pflag),'p'))
        pflag=true;
    else
        pflag=false;
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
    
    % allocate pp if wanted
    if(pflag)
        p=cell(nrecs,max(ncmp));
    else
        depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    end
    
    % detail message
    if(verbose)
        disp('Getting Polynomial Fit to Record(s)');
        print_time_left(0,nrecs);
    end

    % get polynomial and update header if desired
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
            if(pflag)
                for j=1:ncmp(i)
                    p{i,j}=polyfit(time,data(i).dep(:,j),degree(i,j));
                end
            else
                for j=1:ncmp(i)
                    data(i).dep(:,j)=polyval(polyfit(time,...
                        data(i).dep(:,j),degree(i,j)),time);
                end
            end
        % unevenly spaced
        else
            if(pflag)
                for j=1:ncmp(i)
                    p{i,j}=polyfit(double(data(i).ind),data(i).dep(:,j),...
                        degree(i,j));
                end
            else
                for j=1:ncmp(i)
                    data(i).dep(:,j)=polyval(polyfit(...
                        double(data(i).ind),data(i).dep(:,j),...
                        degree(i,j)),double(data(i).ind));
                end
            end
        end
        
        if(~pflag)
            % change class back
            data(i).dep=oclass(data(i).dep);
            
            % adjust header
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    if(pflag)
        data=p;
    else
        % adjust header
        data=changeheader(data,...
            'depmen',depmen,'depmin',depmin,'depmax',depmax);
    end

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
