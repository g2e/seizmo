function [data,roughness]=getspline(data,roughness,ppflag)
%GETSPLINE    Get smooth spline fit of SEIZMO records
%
%    Usage:    data=getspline(data)
%              data=getspline(data,roughness)
%              pp=getspline(data,roughness,'pp')
%              [...,roughness]=getspline(...)
%
%    Description:
%     DATA=GETSPLINE(DATA) returns the smooth spline fit of SEIZMO records
%     in DATA.  For multi-component records, each component is dealt with
%     separately.  To control the roughness of the smoothing spline see the
%     next usage form.
%
%     DATA=GETSPLINE(DATA,ROUGHNESS) sets the roughness of the smoothing
%     spline.  ROUGHNESS should be a value from 0 (linear fit -- the same
%     as REMOVETREND) to 1 (fits all data points thus removing all the data
%     variations).  The default value is determined by the Spline Toolbox
%     function CSAPS and more details about choosing an appropriate value
%     can be found there.  ROUGHNESS can be a scalar or an array of size
%     NRECSxNCMP to specify different values for each record & component.
%
%     PP=GETSPLINE(DATA,ROUGHNESS,'PP') exports the piecewise polynomials
%     for each record & component instead of computing the values.  This
%     could be useful for interpolation.  PP is a NRECSxNCMP cell array of
%     piecewise polynomial structs which may be evaluated using PPVAL.  Any
%     other value for the 3rd argument to GETSPLINE will return the default
%     which is the corresponding values in the SEIZMO struct.
%
%     [...,ROUGHNESS]=GETSPLINE(...) returns the roughness values
%     used for the smooth spline(s).  This is useful for seeing the default
%     values used.  Note that each record and component has its own value.
%
%    Notes:
%     - Requires the Spline Toolbox!
%     - You may also specify using the default roughness using NaN.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Check out the difference various roughness values make:
%     plot0(getspline(data(ones(1,10)),.05:.1:.95));
%
%     % Now do a dry run to get the default roughness
%     % and try a few other roughness values based on that:
%     [sdata,defrough]=getspline(data(1));
%     sdata(2:4)=getspline(data(ones(3,1)),defrough./[100 10 1/10]);
%     plot0(sdata);
%
%     % Return the piecewise polynomial and evaluate it
%     % at 1/10 the sample spacing of the input record:
%     [b,delta,e]=getheader(data(1),'b','delta','e');
%     pp=getspline(data(1),[],'pp');
%     plot1(data(1));
%     hold(gca,'on');
%     plot(b:delta/10:e,ppval(pp{1},b:delta/10:e));
%
%    See also: REMOVESPLINE, GETPOLYNOMIAL, REMOVEMEAN, REMOVETREND,
%              REMOVEPOLYNOMIAL, REMOVEDEADRECORDS, TAPER, CSAPS

%     Version History:
%        Jan. 21, 2014 - initial version
%        Jan. 22, 2014 - ppform output option
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(1,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% default roughness is empty
if(nargin<2 || isempty(roughness)); roughness=nan; end

% default ppflag is empty
if(nargin<3); ppflag=[]; end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt spline fit
try
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % check roughness
    if(~isnumeric(roughness) || ~isreal(roughness) ...
            || any(roughness(:)<0 | roughness>1) ...
            || (~isscalar(roughness) && mod(numel(roughness),nrecs)~=0) ...
            || ~any(size(roughness,1)==[1 nrecs]))
        error('seizmo:getspline:badOrder',...
            'ROUGHNESS must be a scalar or array valued between 0 & 1.');
    end
    if(size(roughness,2)==nrecs && size(roughness,1)==1)
        roughness=roughness(:);
    end
    if(~any(size(roughness,1)==[1 nrecs]))
        error('seizmo:getspline:badOrder',...
            'ROUGHNESS must be a scalar or a NRECSxNCMP array.');
    end
    
    % check ppflag
    if(~isempty(ppflag) && ischar(ppflag) && isequal(lower(ppflag),'pp'))
        ppflag=true;
    else
        ppflag=false;
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
    
    % allocate pp if wanted
    if(ppflag)
        pp=cell(nrecs,max(ncmp));
    else
        depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    end
    
    % detail message
    if(verbose)
        disp('Getting Smooth Spline Fit to Record(s)');
        print_time_left(0,nrecs);
    end

    % get spline and update header if desired
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
            if(ppflag)
                for j=1:ncmp(i)
                    if(isnan(roughness(i,j)))
                        [pp{i,j},roughness(i,j)]=csaps(time,...
                            data(i).dep(:,j));
                    else
                        pp{i,j}=csaps(time,data(i).dep(:,j),...
                            roughness(i,j));
                    end
                end
            else
                for j=1:ncmp(i)
                    if(isnan(roughness(i,j)))
                        [data(i).dep(:,j),roughness(i,j)]=csaps(time,...
                            data(i).dep(:,j),[],time);
                    else
                        data(i).dep(:,j)=csaps(time,data(i).dep(:,j),...
                            roughness(i,j),time);
                    end
                end
            end
        % unevenly spaced
        else
            if(ppflag)
                for j=1:ncmp(i)
                    if(isnan(roughness(i,j)))
                        [pp{i,j},roughness(i,j)]=csaps(...
                            double(data(i).ind),data(i).dep(:,j));
                    else
                        pp{i,j}=csaps(double(data(i).ind),...
                            data(i).dep(:,j),roughness(i,j));
                    end
                end
            else
                for j=1:ncmp(i)
                    if(isnan(roughness(i,j)))
                        [data(i).dep(:,j),roughness(i,j)]=csaps(...
                            double(data(i).ind),data(i).dep(:,j),[],...
                            double(data(i).ind));
                    else
                        data(i).dep(:,j)=csaps(double(data(i).ind),...
                            data(i).dep(:,j),roughness(i,j),...
                            double(data(i).ind));
                    end
                end
            end
        end
        
        if(~ppflag)
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
    
    if(ppflag)
        data=pp;
    else
        % adjust header
        data=changeheader(data,...
            'depmen',depmen,'depmin',depmin,'depmax',depmax);
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
