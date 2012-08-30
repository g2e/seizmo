function [data,dep,ind]=detach(data,option,dpts)
%DETACH    Detach data from SEIZMO records
%
%    Usage:    [data,dep,ind]=detach(data,option,dpts)
%
%    Description:
%     [DATA,DEP,IND]=DETACH(DATA,OPTION,DPTS) detaches DPTS number of
%     points from records in SEIZMO struct DATA and returns them in DEP and
%     IND.  DPTS should be a numeric scalar or an array with 1 element per
%     record.  DEP and IND are cell arrays with each element corresponding
%     to an individual record in DATA.  OPTION must be either 'beginning'
%     or 'ending'.  If OPTION is 'beginning', DEP & IND contain detached
%     data from the beginning of the records and the B & NPTS header fields
%     are updated accordingly.  If OPTION is 'ending', DEP & IND contain
%     detached data from the ending of the records and the E & NPTS header
%     fields are updated accordingly.
%
%    Notes:
%
%    Header changes: NPTS, B, E, DELTA, DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Remove the first 10 points from all records:
%     data=detach(data,'beginning',10);
%
%    See also: ATTACH, CUT

%     Version History:
%        Oct. 10, 2009 - initial version
%        Jan. 26, 2010 - seizmoverbose support
%        Jan. 30, 2010 - minor message update
%        Feb.  2, 2010 - versioninfo caching
%        Mar.  8, 2010 - versioninfo caching dropped
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  2, 2011 - don't allow xyz iftype, doc update, better
%                        checkheader usage, set b/e to undef if detach all
%        Mar. 13, 2012 - seizmocheck fix, use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt detach
try
    % check headers
    data=checkheader(data,'NONTIME_IFTYPE','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % get header info
    [npts,b,delta,e,leven]=getheader(data,...
        'npts','b','delta','e','leven lgc');
    leven=~strcmpi(leven,'false');

    % check option
    validopt={'beginning' 'ending'};
    if(~ischar(option))
        error('seizmo:detach:badOption','OPTION must be a string!');
    end
    if(~strcmpi(option,validopt))
        error('seizmo:detach:badOption',...
            ['OPTION must be one of the following:\n' ...
            sprintf('''%s'' ',validopt{:})]);
    end
    option(1)=upper(option(1));

    % check npts
    if(~isnumeric(dpts) || (~isscalar(dpts) && numel(dpts)~=nrecs))
        error('seizmo:detach:badDPTS',...
            'DPTS must be an integer or an array w/ 1 integer/record!');
    end
    if(isscalar(dpts)); dpts=dpts(ones(nrecs,1),1); end
    if(any(dpts>npts))
        error('seizmo:detach:DPTSgtNPTS',...
            'DPTS must be <= NPTS in each record!')
    end
    
    % detail message
    if(verbose)
        disp(['Detaching Data from Record(s) ' option]);
        print_time_left(0,nrecs);
    end

    % loop over records
    dep=cell(nrecs,1); ind=cell(nrecs,1);
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % unevenly spaced
        if(~leven(i))
            % extract/delete
            switch lower(option)
                case 'beginning'
                    dep{i}=data(i).dep(1:dnpts(i),:);
                    data(i).dep(1:dnpts(i),:)=[];
                    ind{i}=data(i).ind(1:dnpts(i));
                    data(i).ind(1:dnpts(i))=[];
                    npts(i)=npts(i)-dpts(i);
                case 'ending'
                    dep{i}=data(i).dep(end-dnpts(i)+1:end,:);
                    data(i).dep(end-dnpts(i)+1:end,:)=[];
                    ind{i}=data(i).ind(end-dnpts(i)+1:end);
                    data(i).ind(end-dnpts(i)+1:end)=[];
                    npts(i)=npts(i)-dpts(i);
            end

            % timing
            if(npts(i))
                b(i)=data(i).ind(1);
                e(i)=data(i).ind(end);
                if(npts(i)>1); delta=(e(i)-b(i))/(npts(i)-1); end
            else
                b(i)=nan;
                e(i)=nan;
            end
        else % evenly spaced
            % extract/delete
            switch lower(option)
                case 'beginning'
                    dep{i}=data(i).dep(1:dnpts(i),:);
                    data(i).dep(1:dnpts(i),:)=[];
                    npts(i)=npts(i)-dpts(i);
                    b(i)=b(i)+delta(i)*dpts(i);
                case 'ending'
                    dep{i}=data(i).dep(end-dnpts(i)+1:end,:);
                    data(i).dep(end-dnpts(i)+1:end,:)=[];
                    npts(i)=npts(i)-dpts(i);
                    e(i)=e(i)-delta(i)*dpts(i);
            end
            
            % timing
            if(~npts(i))
                b(i)=nan;
                e(i)=nan;
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
