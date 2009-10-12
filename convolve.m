function [data,zf]=convolve(data,tf,zi)
%CONVOLVE    Convolve SEIZMO records with a time function
%
%    Usage:    data=convolve(data,tf)
%              [data,zf]=convolve(data,tf)
%              [data,zf]=convolve(data,tf,zi)
%
%    Description: DATA=CONVOLVE(DATA,TF) convolves a time function TF onto
%     records in SEIZMO struct DATA.  TF must be a numeric vector or a cell
%     array with as many elements as records in DATA with each element
%     being a numeric vector (allows for each record to have a different
%     time function applied).  The sample spacing of the time function is
%     assumed to be the same as that of each record.  Convolution is done
%     in double precision in the time domain using FILTER.  The records are
%     then converted back to their original class.
%
%     [DATA,ZF]=CONVOLVE(DATA,TF) returns the final conditions of the
%     convolution in ZF.  ZF is a cell array with each element
%     corresponding to an individual record in DATA.  So ZF{3} corresponds
%     to DATA(3).  ZF{I} is size [LENGTH(TF)-1 GETNCMP(DATA(I))].  This is
%     useful for convolving individual segments of an extended record
%     without losing information at the segment boundaries.
%
%     [DATA,ZF]=CONVOLVE(DATA,TF,ZI) applies initial conditions ZI to the
%     convolution.  ZI should match the size of the final conditions ZF.
%
%    Notes:
%     - use TIMESHIFT to take care of any delay effects of the convolution
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Convolve a record with an 11-point gaussian:
%      [data1,zf]=convolve(data(1),gausswin(11));
%
%     Now center the gaussian by time shifting by the half-duration:
%      data2=timeshift(data1,-5*getheader(data1,'delta'));
%
%     And plot an overlay to confirm (attaching final conditions too):
%      plot2([data; attach(data2,'append',zf)])
%
%    See also: DECONVOLVE, ATTACH, IIRFILTER, CORRELATE

%     Version History:
%        Oct.  8, 2009 - initial version
%        Oct. 10, 2009 - LEVEN and IFTYPE checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2009 at 20:40 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
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
ncmp=getheader(data);

% cannot do spectral/xyz records
if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
    error('seizmo:convolve:badIFTYPE',...
        'Datatype of records in DATA must be Timeseries or XY!');
end

% cannot do unevenly sampled records
if(any(strcmpi(leven,'false')))
    error('seizmo:convolve:badLEVEN',...
        'Invalid operation on unevenly sampled records!');
end

% number of records
nrecs=numel(data);

% check time function
if(iscell(tf))
    if(numel(tf)~=nrecs)
        error('seizmo:convolve:badTF',...
            'TF must be a cell array with one element per record!');
    end
    ntf=nan(nrecs,1);
    for i=1:nrecs
        if(~isnumeric(tf{i}) || ~isvector(tf{i}) || isempty(tf{i}))
            error('seizmo:convolve:badTF',...
                'TF elements must be non-empty numeric vectors!');
        else
            ntf(i)=numel(tf{i});
            tf{i}=double(tf{i}(:).');
        end
    end
    tfidx=1:nrecs;
else
    if(~isnumeric(tf) || ~isvector(tf) || isempty(tf))
        error('seizmo:convolve:badTF',...
            'TF must be a non-empty numeric vector!');
    end
    ntf=numel(tf);
    tf={double(tf(:).')};
    tfidx=ones(nrecs,1);
end

% check initial conditions
if(nargin>2 && ~isempty(zi))
    if(iscell(zi))
        if(numel(zi)~=nrecs)
            error('seizmo:convolve:badZI',...
                'ZI must be a cell array with one element per record!');
        end
        for i=1:nrecs
            if(~isnumeric(zi{i}) ...
                    || ~isequal(size(zi{i}),[ntf(tfidx(i))-1 ncmp(i)]))
                error('seizmo:convolve:badZI',...
                    ['ZI elements must be numeric arrays ' ...
                    'of size [numel(TF)-1 ncmp]!']);
            else
                zi{i}=double(zi{i});
            end
        end
        ziidx=1:nrecs;
    else
        if(~isnumeric(zi{i}) ...
                || ~isequal(size(zi),[ntf(tfidx(i))-1 ncmp(i)]))
            error('seizmo:convolve:badZI',...
                'ZI must be numeric array of size [numel(TF)-1 ncmp]!');
        end
        zi={double(zi)};
        ziidx=ones(nrecs,1);
    end
end

% loop through records
zf=cell(nrecs,1);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get data class, convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % convolve
    if(nargin==2 || isempty(zi))
        [data(i).dep,zf{i}]=filter(tf{tfidx(i)},1,data(i).dep,[],1);
    else
        [data(i).dep,zf{i}]=...
            filter(tf{tfidx(i)},1,data(i).dep,zi{ziidx(i)},1);
    end
    
    % dep*
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
    
    % restore class
    data(i).dep=oclass(data(i).dep);
    zf{i}=oclass(zf{i});
end

% update header
data=changeheader(data,'depmax',depmax,'depmin',depmin,'depmen',depmen);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
