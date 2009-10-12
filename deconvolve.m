function [data]=deconvolve(data,tf,h2o,zi,zf)
%DECONVOLVE    Spectrally deconvolve a time function from SEIZMO records
%
%    Usage:    data=deconvolve(data,tf)
%              data=deconvolve(data,tf,h2o)
%              data=deconvolve(data,tf,h2o,zi)
%              data=deconvolve(data,tf,h2o,zi,zf)
%
%    Description: DATA=DECONVOLVE(DATA,TF) deconvolves time function TF out
%     of records in SEIZMO struct DATA using spectral division.  TF must be
%     a numeric vector or a cell array with one numeric vector per record
%     in DATA.  The sample spacing of the time function is assumed to be
%     the same as the record it is deconvolved from.  Records in DATA must
%     be evenly spaced and Timeseries or XY datatype.  Note that if TF has
%     spectral amplitudes == 0, then the deconvolution becomes unstable
%     (see the next usage description for a more stable operation).
%
%     DATA=DECONVOLVE(DATA,TF,H2O) adds factor H2O to the spectral
%     amplitudes of TF to avoid division by zero.  This is commonly known
%     as setting the waterlevel.  The best value for H2O varies and often
%     requires some inspection.  The default value is 0 and is the unstable
%     case.  H2O must be a positive scalar or an array of positive values,
%     1 per record.
%
%     DATA=DECONVOLVE(DATA,TF,H2O,ZI) removes the effect of convolution
%     initial conditions ZI from records in DATA before deconvolution.
%     This will remove the effects of the preceeding data that can not be
%     accounted for otherwise.
%
%     DATA=DECONVOLVE(DATA,TF,H2O,ZI,ZF) attaches convolution final
%     conditions ZF to the records in DATA.  This will include the full
%     convolution information for the final points in the record so they
%     are properly accounted for.
%
%    Notes:
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Convolve and then deconvolve a dataset:
%      [data1,zf]=convolve(data,gausswin(13))
%      data2=deconvolve(data1,gausswin(13),0.001,[],zf);
%
%    See also: DECONVOLVETD, CONVOLVE, IIRFILTER, CORRELATE

%     Version History:
%        Oct. 12, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 12, 2009 at 17:50 GMT

% todo:

% check nargin
msg=nargchk(2,5,nargin);
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
    error('seizmo:deconvolve:badIFTYPE',...
        'Datatype of records in DATA must be Timeseries or XY!');
end

% cannot do unevenly sampled records
if(any(strcmpi(leven,'false')))
    error('seizmo:deconvolve:badLEVEN',...
        'Invalid operation on unevenly sampled records!');
end

% number of records
nrecs=numel(data);

% check time function
if(iscell(tf))
    if(numel(tf)~=nrecs)
        error('seizmo:deconvolve:badTF',...
            'TF must be a cell array with one element per record!');
    end
    ntf=nan(nrecs,1);
    for i=1:nrecs
        if(~isnumeric(tf{i}) || ~isvector(tf{i}) || isempty(tf{i}))
            error('seizmo:deconvolve:badTF',...
                'TF elements must be non-empty numeric vectors!');
        elseif(tf{i}(1)==0)
            error('seizmo:deconvolve:badTF',...
                'TF must begin with a non-zero!');
        else
            ntf(i)=numel(tf{i});
            tf{i}=double(tf{i}(:));
        end
    end
    tfidx=1:nrecs;
else
    if(~isnumeric(tf) || ~isvector(tf) || isempty(tf))
        error('seizmo:deconvolve:badTF',...
            'TF must be a non-empty numeric vector!');
    elseif(tf(1)==0)
        error('seizmo:deconvolve:badTF',...
            'TF must begin with a non-zero!');
    end
    ntf=numel(tf);
    tf={double(tf(:))};
    tfidx=ones(nrecs,1);
end

% check water level
if(nargin<3 || isempty(h2o))
    h2o=0;
elseif(~isnumeric(h2o) || (~isscalar(h2o) && numel(h2o)~=nrecs))
    error('seizmo:deconvolve:badH2O',...
        'H2O must be a numeric scalar or an array w/ 1 element/record!');
elseif(any(h2o<0))
    error('seizmo:deconvolve:badH2O','H2O must be positive!');
end
if(isscalar(h2o)); h2o=h2o(ones(nrecs,1),1); end

% check initial conditions
if(nargin>3 && ~isempty(zi))
    if(iscell(zi))
        if(numel(zi)~=nrecs)
            error('seizmo:deconvolve:badZI',...
                'ZI must be a cell array with one element per record!');
        end
        nzi=nan(nrecs,1);
        for i=1:nrecs
            if(~isnumeric(zi{i}) ...
                    || ~isequal(size(zi{i}),[ntf(tfidx(i))-1 ncmp(i)]))
                error('seizmo:deconvolve:badZI',...
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
            error('seizmo:deconvolve:badZI',...
                'ZI must be numeric array of size [numel(TF)-1 ncmp]!');
        end
        nzi=size(zi,1);
        zi={double(zi)};
        ziidx=ones(nrecs,1);
    end
end

% check initial conditions
if(nargin>4 && ~isempty(zf))
    if(iscell(zf))
        if(numel(zf)~=nrecs)
            error('seizmo:deconvolve:badZF',...
                'ZF must be a cell array with one element per record!');
        end
        for i=1:nrecs
            if(~isnumeric(zf{i}) ...
                    || ~isequal(size(zf{i}),[ntf(tfidx(i))-1 ncmp(i)]))
                error('seizmo:deconvolve:badZF',...
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
            error('seizmo:deconvolve:badZF',...
                'ZF must be numeric array of size [numel(TF)-1 ncmp]!');
        end
        zf={double(zf)};
        zfidx=ones(nrecs,1);
    end
end

% loop through records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get data class, convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % clear initial conditions
    if(nargin>3 && ~isempty(zi))
        data(i).dep(1:nzi(ziidx(i)),:)=...
            data(i).dep(1:nzi(ziidx(i)),:)-zi{ziidx(i)};
    end
    
    % waterlevel stabilized spectral division
    if(nargin>3 && ~isempty(zi))
        nspts=2^(nextpow2(npts(i)+ntf(tfidx(i))-1)+1);
        tmp=fft(tf{tfidx(i)}(:,ones(1,ncmp(i))),nspts);
        tmp=(abs(tmp)+h2o(i)).*exp(j*angle(tmp));
        tmp=ifft(fft([data(i).dep; zf{zfidx(i)}],nspts)./tmp,'symmetric');
        data(i).dep=tmp(1:npts(i),:);
    else
        nspts=2^(nextpow2(npts(i))+1);
        tmp=fft(tf{tfidx(i)}(:,ones(1,ncmp(i))),nspts);
        tmp=(abs(tmp)+h2o(i)).*exp(j*angle(tmp));
        tmp=ifft(fft(data(i).dep,nspts)./tmp,'symmetric');
        data(i).dep=tmp(1:npts(i),:);
    end
    
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
