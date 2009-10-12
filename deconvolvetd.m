function [data,residue]=deconvolvetd(data,tf,zi,zf)
%DECONVOLVETD    Time-domain Deconvolution of SEIZMO records
%
%    Usage:    data=deconvolvetd(data,tf)
%              [data,residue]=deconvolvetd(data,tf)
%              [...]=deconvolvetd(data,tf,zi)
%              [...]=deconvolvetd(data,tf,zi,zf)
%
%    Description: DATA=DECONVOLVETD(DATA,TF) deconvolves impulse response
%     TF out of records in SEIZMO struct DATA.  TF must be a numeric vector
%     or a cell array with one numeric vector per record in DATA.  The
%     sample spacing of the time function is assumed to be the same as the
%     record it is deconvolved from.  Note that deconvolution is done in
%     double precision but the operation is highly unstable (no waterlevel
%     factor or something of that sort is used).  In general, this function
%     is only useful for the simplest operations.
%
%     [DATA,RESIDUE]=DECONVOLVETD(DATA,TF) returns the residue from the
%     deconvolution as RESIDUE.  RESIDUE is a cell array with a numeric
%     array per record in DATA.  The inputs and outputs are thus basically
%     related by DATAIN = convolve(DATAOUT,TF) + RESIDUE.
%
%     [...]=DECONVOLVETD(DATA,TF,ZI) removes the effect of convolution
%     initial conditions ZI from records in DATA before deconvolution.
%     This will remove the effects of the preceeding data that can not be
%     accounted for otherwise.
%
%     [...]=DECONVOLVETD(DATA,TF,ZI,ZF) attaches convolution final
%     conditions ZF to the records in DATA.  This will include the full
%     convolution information for the final points in the record so they
%     are properly accounted for.
%
%    Notes:
%     - The algorithm in DECONVOLVETD is unstable and is not recommended.
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Convolve a gaussian onto a record of delta functions and then attempt
%     to deconvolve the gaussian from the convolved record:
%      data0=bseizmo(1:100,[0 0 0 1 0 0 0 0 -2 0 0 0 0 1.5 zeros(1,86)]);
%      data1=convolve(data0,gausswin(9));
%      data2=deconvolvetd(data1,gausswin(9));
%      plot2([data0 data2],'xlimits',[-4 50]); % dont plot bad part
%
%    See also: DECONVOLVE, CONVOLVE, IIRFILTER, CORRELATE

%     Version History:
%        Oct.  9, 2009 - initial version
%        Oct. 10, 2009 - LEVEN and IFTYPE checks
%        Oct. 12, 2009 - name changed to DECONVOLVETD
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 12, 2009 at 06:55 GMT

% todo:

% check nargin
msg=nargchk(2,4,nargin);
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
[npts,ncmp]=getheader(data,'npts','ncmp');

% cannot do spectral/xyz records
if(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
    error('seizmo:deconvolvetd:badIFTYPE',...
        'Datatype of records in DATA must be Timeseries or XY!');
end

% cannot do unevenly sampled records
if(any(strcmpi(leven,'false')))
    error('seizmo:deconvolvetd:badLEVEN',...
        'Invalid operation on unevenly sampled records!');
end

% number of records
nrecs=numel(data);

% check time function
if(iscell(tf))
    if(numel(tf)~=nrecs)
        error('seizmo:deconvolvetd:badTF',...
            'TF must be a cell array with one element per record!');
    end
    ntf=nan(nrecs,1);
    for i=1:nrecs
        if(~isnumeric(tf{i}) || ~isvector(tf{i}) || isempty(tf{i}))
            error('seizmo:deconvolvetd:badTF',...
                'TF elements must be non-empty numeric vectors!');
        elseif(tf{i}(1)==0)
            error('seizmo:deconvolvetd:badTF',...
                'TF must begin with a non-zero!');
        else
            ntf(i)=numel(tf{i});
            tf{i}=double(tf{i}(:).');
        end
    end
    tfidx=1:nrecs;
else
    if(~isnumeric(tf) || ~isvector(tf) || isempty(tf))
        error('seizmo:deconvolvetd:badTF',...
            'TF must be a non-empty numeric vector!');
    elseif(tf(1)==0)
        error('seizmo:deconvolvetd:badTF',...
            'TF must begin with a non-zero!');
    end
    ntf=numel(tf);
    tf={double(tf(:).')};
    tfidx=ones(nrecs,1);
end

% check initial conditions
if(nargin>2 && ~isempty(zi))
    if(iscell(zi))
        if(numel(zi)~=nrecs)
            error('seizmo:deconvolvetd:badZI',...
                'ZI must be a cell array with one element per record!');
        end
        nzi=nan(nrecs,1);
        for i=1:nrecs
            if(~isnumeric(zi{i}) ...
                    || ~isequal(size(zi{i}),[ntf(tfidx(i))-1 ncmp(i)]))
                error('seizmo:deconvolvetd:badZI',...
                    ['ZI elements must be numeric arrays ' ...
                    'of size [numel(TF)-1 ncmp]!']);
            else
                nzi(i)=size(zi{i},1);
                zi{i}=double(zi{i});
            end
        end
        ziidx=1:nrecs;
    else
        if(~isnumeric(zi{i}) ...
                || ~isequal(size(zi),[ntf(tfidx(i))-1 ncmp(i)]))
            error('seizmo:deconvolvetd:badZI',...
                'ZI must be numeric array of size [numel(TF)-1 ncmp]!');
        end
        nzi=size(zi,1);
        zi={double(zi)};
        ziidx=ones(nrecs,1);
    end
end

% check initial conditions
if(nargin>3 && ~isempty(zf))
    if(iscell(zf))
        if(numel(zf)~=nrecs)
            error('seizmo:deconvolvetd:badZF',...
                'ZF must be a cell array with one element per record!');
        end
        for i=1:nrecs
            if(~isnumeric(zf{i}) ...
                    || ~isequal(size(zf{i}),[ntf(tfidx(i))-1 ncmp(i)]))
                error('seizmo:deconvolvetd:badZF',...
                    ['ZF elements must be numeric arrays ' ...
                    'of size [numel(TF)-1 ncmp]!']);
            else
                zf{i}=double(zf{i});
            end
        end
        zfidx=1:nrecs;
    else
        if(~isnumeric(zf{i}) ...
                || ~isequal(size(zf),[ntf(tfidx(i))-1 ncmp(i)]))
            error('seizmo:deconvolvetd:badZF',...
                'ZF must be numeric array of size [numel(TF)-1 ncmp]!');
        end
        zf={double(zf)};
        zfidx=ones(nrecs,1);
    end
end

% loop through records
residue=cell(nrecs,1);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get data class, convert to double precision
    oclassstr=class(data(i).dep);
    oclass=str2func(oclassstr);
    data(i).dep=double(data(i).dep);
    
    % clear initial conditions
    if(nargin>2 && ~isempty(zi))
        data(i).dep(1:nzi(ziidx(i)),:)=...
            data(i).dep(1:nzi(ziidx(i)),:)-zi{ziidx(i)};
    end
    
    % deconvolve
    temp=cell(1,ncmp(i));
    if(nargin>3 && ~isempty(zf))
        for j=1:ncmp(i)
            [data(i).dep(:,j),temp{j}]=filter(...
                [data(i).dep(:,j); zf{zfidx(i)}(:,j)],...
                tf{tfidx(i)},[1; zeros(npts(i)-1,1)]);
        end
    else
        for j=1:ncmp(i)
            % adds in zeros as final conditions
            [data(i).dep(:,j),temp{j}]=filter(...
                [data(i).dep(:,j); zeros(ntf(tfidx(i))-1,1)],...
                tf{tfidx(i)},[1; zeros(npts(i)-1,1)]);
        end
    end
    
    % get residue
    temp=cell2mat(temp);
    residue{i}=zeros(npts(i)+ntf(tfidx(i))-1,ncmp(i),oclassstr);
    residue{i}(npts(i)+1:end,:)=tf{tfidx(i)}(1).*temp(1:ntf(tfidx(i))-1,:);
    
    % dep*
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
    
    % restore class
    data(i).dep=oclass(data(i).dep);
end

% update header
data=changeheader(data,'depmax',depmax,'depmin',depmin,'depmen',depmen);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
