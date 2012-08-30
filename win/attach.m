function [data]=attach(data,option,dep,ind)
%ATTACH    Attach data to SEIZMO records
%
%    Usage:    data=attach(data,option,dep)
%              data=attach(data,option,dep,ind)
%
%    Description:
%     DATA=ATTACH(DATA,OPTION,DEP) attaches DEP to the dependent data in
%     records in SEIZMO struct DATA.  DEP must be a numeric array or a cell
%     array of numeric arrays (one per record).  Numeric arrays in DEP must
%     have the same number of columns as the number of components in the
%     record it is being attached to.  OPTION is either 'beginning' or
%     'ending' and decides how DEP is attached to the dependent data.  For
%     instance, if OPTION is 'ending', DEP is attached at the end and the E
%     header field is adjusted.  If OPTION is 'beginning', DEP is attached
%     at the beginning of the dependent dataset and the B header field is
%     adjusted accordingly.
%
%     DATA=ATTACH(DATA,OPTION,DEP,IND) attaches IND to the independent
%     component data.  IND must be supplied if there is any unevenly
%     sampled records.  IND is ignored for evenly sampled records.  The
%     number of rows in IND and DEP numeric arrays must match.
%
%    Notes:
%
%    Header changes: NPTS, B, E, DELTA, DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Attach final conditions to convolved records:
%     [data,zf]=convolve(data,gausswin(11));
%     data=attach(data,'ending',zf);
%
%     % Add in points to unevenly sampled record (later sorting by timing):
%     data=checkheader(...
%         attach(data,'ending',dep,ind),'nonmonotonic_ind','fix')
%
%    See also: DETACH, CUT

%     Version History:
%        Oct. 10, 2009 - initial version
%        Oct. 13, 2009 - minor doc update
%        Oct. 26, 2009 - added 'prepend' and 'append' as valid options,
%                        fixed bug in matrix entry
%        Jan. 26, 2010 - seizmoverbose support
%        Feb.  2, 2010 - versioninfo caching
%        Mar.  8, 2010 - versioninfo caching dropped
%        Aug. 16, 2010 - nargchk fix, better checkheader utilization
%        Nov.  2, 2011 - doc update
%        Mar. 13, 2012 - seizmocheck fix, use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 12:45 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt attach
try
    % check headers
    data=checkheader(data,'NONTIME_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % get header info
    [npts,ncmp,b,delta,e,leven]=getheader(data,...
        'npts','ncmp','b','delta','e','leven lgc');
    leven=~strcmpi(leven,'false');

    % check option
    validopt={'beginning' 'ending' 'prepend' 'append'};
    if(~ischar(option))
        error('seizmo:attach:badOption','OPTION must be a string!');
    end
    if(~strcmpi(option,validopt))
        error('seizmo:attach:badOption',...
            ['OPTION must be one of the following:\n' ...
            sprintf('''%s'' ',validopt{:})]);
    end
    option(1)=upper(option(1));

    % check dep
    if(iscell(dep))
        if(numel(dep)~=nrecs)
            error('seizmo:attach:badDEP',...
                'DEP must be a cell array with one element per record!');
        end
        ndep=nan(nrecs,1);
        for i=1:nrecs
            if(~isnumeric(dep{i}) ...
                    || (~isempty(dep{i}) && size(dep{i},2)~=ncmp(i)))
                error('seizmo:attach:badDEP',...
                    'DEP elements must be real arrays of size [n ncmp]!');
            end
            ndep(i)=size(dep{i},1);
        end
        depidx=1:nrecs;
    else
        if(~isnumeric(dep) || (~isempty(dep) && any(size(dep,2)~=ncmp)))
            error('seizmo:attach:badDEP',...
                'DEP must be a real array of size [n ncmp]!');
        end
        ndep=size(dep,1);
        depidx=ones(nrecs,1);
        dep={dep};
    end

    % check ind
    if(nargin==3 && any(~leven))
        error('seizmo:attach:INDnecessary',...
            ['IND argument must be given if ' ...
            'DATA has unevenly sampled records!']);
    elseif(nargin>3 && ~isempty(ind))
        if(iscell(ind))
            if(numel(ind)~=nrecs)
                error('seizmo:attach:badIND',...
                    'IND must be a cell array with 1 element per record!');
            end
            nind=nan(nrecs,1);
            for i=1:nrecs
                if(~isnumeric(ind{i}) ...
                        || (~isempty(ind{i}) && size(ind{i},2)~=1))
                    error('seizmo:attach:badIND',...
                        ['IND elements must be numeric ' ...
                        'arrays of size [n 1]!']);
                end
                nind(i)=size(ind{i},1);
            end
            indidx=1:nrecs;
        else
            if(~isnumeric(ind) || any(size(ind,2)~=1))
                error('seizmo:attach:badIND',...
                    'IND must be a numeric array of size [n 1]!');
            end
            nind=size(ind,1);
            indidx=ones(nrecs,1);
            ind={ind};
        end
        if(~isequal(nind(indidx(1:nrecs)),ndep(depidx(1:nrecs))))
            error('seizmo:attach:IndDepSizeMismatch',...
                'DEP & IND must have the same number of rows!');
        end
    end
    
    % detail message
    if(verbose)
        disp(['Attaching Data to Record(s) ' option]);
        print_time_left(0,nrecs);
    end

    % loop over records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % get data class
        oclass=str2func(class(data(i).dep));

        % unevenly spaced
        if(~leven(i))
            % reclass attachment
            tmpd=oclass(dep{depidx(i)});
            oclassi=str2func(class(data(i).ind));
            tmpi=oclassi(ind{indidx(i)});

            % attach
            switch lower(option)
                case {'beginning' 'prepend'}
                    data(i).dep=[tmpd; data(i).dep];
                    data(i).ind=[tmpi; data(i).ind];
                case {'ending' 'append'}
                    data(i).dep=[data(i).dep; tmpd];
                    data(i).ind=[data(i).ind; tmpi];
            end

            % npts/timing
            npts(i)=npts(i)+ndep(depidx(i));
            if(npts(i))
                b(i)=data(i).ind(1);
                e(i)=data(i).ind(end);
                if(npts(i)>1); delta=(e(i)-b(i))/(npts(i)-1); end
            else
                b(i)=nan;
                e(i)=nan;
            end
        else % evenly spaced
            % reclass attachment
            tmpd=oclass(dep{depidx(i)});

            % npts
            npts(i)=npts(i)+ndep(depidx(i));

            % attach/timing
            switch lower(option)
                case {'beginning' 'prepend'}
                    data(i).dep=[tmpd; data(i).dep];
                    b(i)=b(i)-delta(i)*ndep(depidx(i));
                case {'ending' 'append'}
                    data(i).dep=[data(i).dep; tmpd];
                    e(i)=e(i)+delta(i)*ndep(depidx(i));
            end
        end

        % dep*
        if(npts(i))
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,'npts',npts,'b',b,'delta',delta,'e',e,...
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
