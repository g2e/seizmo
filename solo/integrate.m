function [data]=integrate(data,option)
%INTEGRATE    Integrate SEIZMO records
%
%    Usage:    data=integrate(data)
%              data=integrate(data,method)
%
%    Description:
%     INTEGRATE(DATA) will integrate SEIZMO records in DATA using the
%     assumption that the record's points are at the midpoints of the
%     integrated record and that their values give the difference between
%     those points.  This operation will add 1 point to the records,
%     starting with the first point set to 0.  This also is attempted for
%     uneven records but may fail if the midpoint assumption cannot be true
%     (the independent component must increase monotonically).  The
%     trapezoidal rule is used as a fallback in this instance (a warning is
%     issued as well).
%
%     INTEGRATE(DATA,METHOD) allows choosing the integration method.  The
%     methods allowed are:
%       'rectangular' - rectangular rule, starts with zero
%       'trapezoidal' - trapezoidal rule, starts with zero
%       'rectangular-sac' - SAC compatible
%       'trapezoidal-sac' - SAC compatible
%       'midpoint' - (default) undoes DIFFERENTIATE - offset rectangular
%
%    Notes:
%     - rectangular adds 1 point if the record is even, none if not and the
%       timing is left the same (e increases by delta for even records)
%     - trapezoidal does not change the timing of the record
%     - rectangular-sac does not change the timing of the record
%     - trapezoidal-sac adjusts the timing to the midpoints, losing a point
%     - midpoint adds 1 point and shifts times to the midpoints (first
%       point one half sample interval earlier - this is nontrivial for
%       an unevenly spaced records but a good estimation is made with the
%       the midpoint assumption)
%
%    Header changes: B, E, NPTS, DELTA, DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     % Check how good integrate undoes differentiate:
%     plot1(subtractrecords(data,integrate(differentiate(data))))
%
%    See also: DIFFERENTIATE, DIVIDEOMEGA, MULTIPLYOMEGA

%     Version History:
%        Nov. 12, 2008 - initial version
%        Nov. 25, 2008 - finally combined integrt, integrt2 and SAC int
%                        methods into one function, fixed several bugs too
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        May   8, 2009 - uses expanded idep unit set
%        Sep.  8, 2009 - drop SWAP for DEAL
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging, forced dim on some functions
%        May   6, 2010 - fixed ncmp bug, slimmer code for units exchange,
%                        fixed interaction of midpoint warning & scrollbar
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 13, 2012 - doc update, seizmocheck fix, better checkheader
%                        usage, use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt integration
try
    % check data header
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check option
    if(nargin==1 || isempty(option))
        option='midpoint';
    elseif(~ischar(option) && ~iscellstr(option))
        error('seizmo:integrate:badOption',...
            'OPTION must be a char or cellstr array!');
    end

    % expand scalar option
    option=cellstr(option);
    nopt=numel(option);
    if(nopt==1)
        option(1:nrecs,1)=option;
    elseif(nopt~=nrecs)
        error('seizmo:integrate:badOption',...
            'OPTION must be a single option or one option per record!');
    end

    % get header info
    [b,e,delta,npts,ncmp,leven,idep]=getheader(data,...
        'b','e','delta','npts','ncmp','leven lgc','idep id');
    leven=~strcmpi(leven,'false');
    
    % detail message
    if(verbose)
        disp('Integrating Record(s)');
        print_time_left(0,nrecs);
    end

    % integrate
    [depmen,depmin,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % skip dataless
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % which integration style?
        switch option{i}
            % slope at point gives offset of next point from the point
            %
            % start with zero at the first point, and shoot to npts+1
            % or just to npts for uneven (we can't get the time of npts+1)
            case {'rect' 'rectangular'}
                if (leven(i))
                    data(i).dep=[zeros(1,ncmp(i)); ...
                        delta(i)*cumsum(data(i).dep,1)];
                    e(i)=e(i)+delta(i);
                    npts(i)=npts(i)+1;
                else
                    dt=diff(double(data(i).ind),1,1);
                    data(i).dep(1:end-1,:)=data(i).dep(1:end-1,:)...
                        .*dt(:,ones(1,ncmp(i)));
                    data(i).dep=[zeros(1,ncmp(i)); cumsum(data(i).dep,1)];
                    data(i).dep(end,:)=[];
                end
            case {'rect-sac' 'rectangular-sac'}
                % slope at point gives offset of point from the last point
                %
                % don't add 1st point (it's zero) to keep npts the same
                if (leven(i))
                    data(i).dep=delta(i)*cumsum(data(i).dep,1);
                else
                    dt=diff(double(data(i).ind),1,1);
                    data(i).dep(2:end,:)=data(i).dep(2:end,:)...
                        .*dt(:,ones(1,ncmp(i)));
                    data(i).dep=cumsum(data(i).dep,1);
                end
            case {'trap' 'trapezoidal'}
                % area of the trapezoid between input points i and i+1
                % gives the offset of output point i+1 from i
                %
                % start first output point at time of first input point
                % setting it to zero, as we have to start somewhere
                if (leven(i))
                    data(i).dep=delta(i)*cumtrapz(data(i).dep,1);
                else
                    data(i).dep=...
                        cumtrapz(double(data(i).ind),data(i).dep,1);
                end
            case {'trap-sac' 'trapezoidal-sac'}
                % SAC drops the first point (is zero)
                % and shifts times to midpoints
                %
                % this is a -delta/2 time shift from trapezoidal
                if (leven(i))
                    data(i).dep=delta(i)*cumtrapz(data(i).dep,1);
                    data(i).dep(1,:)=[];
                    b(i)=b(i)+delta(i)/2;
                    e(i)=e(i)-delta(i)/2;
                    npts(i)=npts(i)-1;
                else
                    data(i).ind=double(data(i).ind);
                    data(i).dep=cumtrapz(data(i).ind,data(i).dep,1);
                    data(i).dep(1,:)=[];
                    data(i).ind(2:end)=...
                        data(i).ind(2:end)-diff(data(i).ind(:),1,1)/2;
                    data(i).ind(1)=[];
                    b(i)=data(i).ind(1);
                    e(i)=data(i).ind(end);
                    data(i).ind=oclass(data(i).ind);
                    npts(i)=npts(i)-1;
                    delta(i)=(e(i)-b(i))/(npts(i)-1);
                end
            case 'midpoint'
                % input points represent slopes at midpoints
                % of output series
                %
                % setting the first output point (half delta
                % before the first input point) to zero, we
                % use the slope to get to the subsequent points
                if(leven(i))
                    data(i).dep=[zeros(1,ncmp(i)); ...
                        delta(i)*cumsum(data(i).dep,1)];
                    b(i)=b(i)-delta(i)/2;
                    e(i)=e(i)+delta(i)/2;
                    npts(i)=npts(i)+1;
                else
                    % Assume each time point represents a midpoint between
                    % the original time points - get the original times.
                    %
                    % 1 more unknown than knowns
                    % ie. (n midpoints, n+1 original points)
                    %
                    % Solve for available solution range under condition
                    % that all time intervals are positive.  Midpoint of
                    % range gives a good solution.
                    %
                    % Range is just smallest time interval.  This is the
                    % lack of knowledge in this case.
                    data(i).ind=double(data(i).ind);
                    dt=diff(data(i).ind,1,1);
                    [range,pos]=min(dt,[],1);

                    % quick check that the midpoint assumption is plausible
                    if(range<=0 ...
                            || all((dt(2:end-1)-dt(1:end-2))>dt(3:end)) ...
                            || all(dt(1:end-2)<(dt(2:end-1)-dt(3:end))));
                        % nope - trapezoidal rule
                        warning('seizmo:integrate:failedAssumption',...
                            ['Midpoint assumption failed for record %d.'...
                            '\nUsing the trapezoidal rule instead.'],i);
                        data(i).dep=...
                            oclass(cumtrapz(data(i).ind,data(i).dep,1));
                        data(i).ind=oclass(data(i).ind);
                        depmen(i)=nanmean(data(i).dep(:));
                        depmin(i)=min(data(i).dep(:));
                        depmax(i)=max(data(i).dep(:));
                        % detail message
                        if(verbose); print_time_left(i,nrecs,true); end
                        continue;
                    end

                    % find original times
                    orig=zeros(npts(i)+1,1);
                    orig(pos+1)=data(i).ind(pos)+range/2;
                    for j=pos+1:npts(i)
                        orig(j+1)=2*data(i).ind(j)-orig(j);
                    end
                    for j=pos:-1:1
                        orig(j)=2*data(i).ind(j)-orig(j+1);
                    end

                    % check monotonicity of times
                    dto=diff(orig);
                    if(min(dto)<=0)
                        % non-monotonic - bail out, use trapezoidal rule
                        warning('seizmo:integrate:failedAssumption',...
                            ['Midpoint assumption failed for record %d.'...
                            '\nUsing the trapezoidal rule instead.'],i);
                        data(i).dep=...
                            oclass(cumtrapz(data(i).ind,data(i).dep,1));
                        data(i).ind=oclass(data(i).ind);
                        depmen(i)=nanmean(data(i).dep(:));
                        depmin(i)=min(data(i).dep(:));
                        depmax(i)=max(data(i).dep(:));
                        % detail message
                        if(verbose); print_time_left(i,nrecs,true); end
                        continue;
                    end

                    % cumsum
                    data(i).ind=oclass(orig);
                    data(i).dep=[zeros(1,ncmp(i));...
                        cumsum(dto(:,ones(1,ncmp(i))).*data(i).dep,1)];
                    npts(i)=npts(i)+1;
                    b(i)=data(i).ind(1);
                    e(i)=data(i).ind(end);
                    delta(i)=(e(i)-b(i))/(npts(i)-1);
                end
            otherwise
                error('seizmo:integrate:badOption','Unknown OPTION!');
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % get dependent component values
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % change dependent component type
    newunit={'icrackle' 'isnap' 'ijerk' 'iacc' 'ivel' 'idisp' 'iabsmnt' ...
        'iabsity' 'iabseler' 'iabserk' 'iabsnap' 'iabsackle' 'iabspop'};
    [tf,idx]=ismember(idep,[{'ipop'} newunit(1:end-1)]);
    idep(tf)=newunit(idx(tf));
    idep(~tf)={'iunkn'};

    % update headers
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax,...
        'delta',delta,'b',b,'e',e,'npts',npts,'idep',idep);

    % toggle struct checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
