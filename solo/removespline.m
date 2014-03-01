function [data]=removespline(data,roughness)
%REMOVESPLINE    Remove smooth spline fit from SEIZMO records
%
%    Usage:    data=removespline(data)
%              data=removespline(data,roughness)
%              [data,roughness]=removespline(data,...)
%
%    Description:
%     DATA=REMOVESPLINE(DATA) removes a smooth spline fit from SEIZMO
%     records.  For multi-component records, each component is dealt with
%     separately.  To control the roughness of the smoothing spline see the
%     next usage form.
%
%     DATA=REMOVESPLINE(DATA,ROUGHNESS) sets the roughness of the smoothing
%     spline.  ROUGHNESS should be a value from 0 (linear fit -- the same
%     as REMOVETREND) to 1 (fits all data points thus removing all the data
%     variations).  The default value is determined by the Spline Toolbox
%     function CSAPS and more details about choosing an appropriate value
%     can be found there.  ROUGHNESS can be a scalar or an array of size
%     NRECSxNCMP to specify different values for each record & component.
%
%     [DATA,ROUGHNESS]=REMOVESPLINE(DATA,...) returns the roughness values
%     used for the smooth spline removal.  This is useful for seeing the
%     default values used.  Note that each record and component has its own
%     value.
%
%    Notes:
%     - Requires the Spline Toolbox!
%     - You may also specify using the default roughness using NaN.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Check out the difference various roughness values make:
%     plot1(removespline(data(ones(1,10)),.9:.01:.99))
%
%    See also: GETSPLINE, REMOVEPOLYNOMIAL, REMOVEMEAN, REMOVETREND,
%              GETPOLYNOMIAL, REMOVEDEADRECORDS, TAPER, CSAPS

%     Version History:
%        Jan. 21, 2014 - initial version
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% default roughness is empty
if(nargin<2); roughness=nan; end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt spline removal
try
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check roughness
    if(~isnumeric(roughness) || ~isreal(roughness) ...
            || any(roughness(:)<0 | roughness(:)>1) ...
            || (~isscalar(roughness) && mod(numel(roughness),nrecs)~=0) ...
            || ~any(size(roughness,1)==[1 nrecs]))
        error('seizmo:removespline:badOrder',...
            'ROUGHNESS must be a scalar or array valued between 0 & 1.');
    end
    if(size(roughness,2)==nrecs && size(roughness,1)==1)
        roughness=roughness(:);
    end
    if(~any(size(roughness,1)==[1 nrecs]))
        error('seizmo:removespline:badOrder',...
            'ROUGHNESS must be a scalar or a NRECSxNCMP array.');
    end

    % header info
    [delta,npts,ncmp,leven]=getheader(data,...
        'delta','npts','ncmp','leven lgc');
    leven=~strcmpi(leven,'false');
    
    % expand roughness
    if(isscalar(roughness))
        roughness=roughness(ones(nrecs,max(ncmp)));
    elseif(size(roughness,2)==1 && max(ncmp)>1)
        roughness=roughness(:,ones(1,max(ncmp)));
    elseif(size(roughness,2)==max(ncmp) && size(roughness,1)==1)
        roughness=roughness(ones(nrecs,1),:);
    end
    
    % detail message
    if(verbose)
        disp('Removing Smooth Spline Fit from Record(s)');
        print_time_left(0,nrecs);
    end

    % remove spline and update header
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
                if(isnan(roughness(i,j)))
                    [tmp,roughness(i,j)]=csaps(time,data(i).dep(:,j),[],...
                        time);
                    data(i).dep(:,j)=data(i).dep(:,j)-tmp;
                else
                    data(i).dep(:,j)=data(i).dep(:,j)-csaps(time,...
                        data(i).dep(:,j),roughness(i,j),time);
                end
            end
        % unevenly spaced
        else
            for j=1:ncmp(i)
                if(isnan(roughness(i,j)))
                    [tmp,roughness(i,j)]=csaps(double(data(i).ind),...
                        data(i).dep(:,j),[],double(data(i).ind));
                    data(i).dep(:,j)=data(i).dep(:,j)-tmp;
                else
                    data(i).dep(:,j)=data(i).dep(:,j)-csaps(...
                        double(data(i).ind),data(i).dep(:,j),...
                        roughness(i,j),double(data(i).ind));
                end
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
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
