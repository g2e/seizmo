function [data]=attach(data,option,dep,ind)
%ATTACH    Attach data to SEIZMO records
%
%    Usage:    data=attach(data,option,dep)
%              data=attach(data,option,dep,ind)
%
%    Description: DATA=ATTACH(DATA,OPTION,DEP) attaches DEP to the
%     dependent data in records in SEIZMO struct DATA.  DEP must be a
%     numeric array or a cell array of numeric arrays, one per record in
%     DATA.  Numeric arrays in DEP must have the same number of columns as
%     the number of components in the record it is being attached to.
%     OPTION is either 'beginning' or 'ending' and decides how DEP is
%     attached to the dependent data.  For instance, if OPTION is 'ending'
%     then DEP is attached at the end and the E header field is adjusted.
%     If OPTON is 'beginning', DEP is attached at the beginning of the
%     dependent dataset and the B header field is adjusted accordingly.
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
%     Attach final conditions to convolved records:
%      [data,zf]=convolve(data,gausswin(11));
%      data=attach(data,'ending',zf);
%
%     Add in points to unevenly sampled record (later sorting by timing):
%      data=checkheader(...
%          attach(data,'ending',dep,ind),'nonmonotonic_ind','fix')
%
%    See also: DETACH, CUT

%     Version History:
%        Oct. 10, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2009 at 21:40 GMT

% todo:

% check nargin
msg=nargchk(3,4,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% get header info
leven=getlgc(data,'leven');
iftype=getenumid(data,'iftype');
[npts,ncmp,b,delta,e]=getheader(data,'npts','ncmp','b','delta','e');

% cannot do spectral/xyz records
if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
    error('seizmo:attach:badIFTYPE',...
        'Datatype of records in DATA must be Timeseries or XY!');
end

% number of records
nrecs=numel(data);

% check option
validopt={'beginning' 'ending'};
if(~ischar(option))
    error('seizmo:attach:badOption','OPTION must be a string!');
end
if(~strcmpi(option,validopt))
    error('seizmo:attach:badOption',...
        ['OPTION must be one of the following:\n' ...
        sprintf('''%s'' ',validopt{:})]);
end

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
                'DEP elements must be numeric arrays of size [n ncmp]!');
        end
        ndep(i)=size(dep{i},1);
    end
    depidx=1:nrecs;
else
    if(~isnumeric(dep) || (~isempty(dep) && any(size(dep,2)~=ncmp)))
        error('seizmo:attach:badDEP',...
            'DEP must be a numeric array of size [n ncmp]!');
    end
    dep={dep};
    ndep=size(dep,1);
    depidx=ones(nrecs,1);
end

% check ind
if(nargin==3 && any(strcmpi(leven,'false')))
    error('seizmo:attach:INDnecessary',...
        ['IND argument must be given if ' ...
        'DATA has unevenly sampled records!']);
elseif(nargin>3 && ~isempty(ind))
    if(iscell(ind))
        if(numel(ind)~=nrecs)
            error('seizmo:attach:badIND',...
                'IND must be a cell array with one element per record!');
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
        ind={ind};
        nind=size(ind,1);
        indidx=ones(nrecs,1);
    end
    if(~isequal(nind(indidx(1:nrecs)),ndep(depidx(1:nrecs))))
        error('seizmo:attach:IndDepSizeMismatch',...
            'DEP & IND must have the same number of rows!');
    end
end

% loop over records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % get data class
    oclass=str2func(class(data(i).dep));
    
    % unevenly spaced
    if(strcmpi(leven{i},'false'))
        % reclass attachment
        tmpd=oclass(dep{depidx(i)});
        oclassi=str2func(class(data(i).ind));
        tmpi=oclassi(ind{indidx(i)});
        
        % attach
        switch lower(option)
            case 'beginning'
                data(i).dep=[tmpd; data(i).dep];
                data(i).ind=[tmpi; data(i).ind];
            case 'ending'
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
            case 'beginning'
                data(i).dep=[tmpd; data(i).dep];
                b(i)=b(i)-delta(i)*ndep(depidx(i));
            case 'ending'
                data(i).dep=[data(i).dep; tmpd];
                e(i)=e(i)+delta(i)*ndep(depidx(i));
        end
    end
    
    % dep*
    if(npts(i))
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    end
end

% update header
data=changeheader(data,'npts',npts,'b',b,'delta',delta,'e',e,...
    'depmax',depmax,'depmin',depmin,'depmen',depmen);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
