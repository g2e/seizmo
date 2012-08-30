function [data]=deconvolve(data,tf,delay,h2o,frange,zi,zf)
%DECONVOLVE    Spectrally deconvolve a time function from SEIZMO records
%
%    Usage:    data=deconvolve(data,tf)
%              data=deconvolve(data,tf,delay)
%              data=deconvolve(data,tf,delay,h2o)
%              data=deconvolve(data,tf,delay,h2o,frange)
%              data=deconvolve(data,tf,delay,h2o,frange,zi)
%              data=deconvolve(data,tf,delay,h2o,frange,zi,zf)
%
%    Description:
%     DATA=DECONVOLVE(DATA,TF) deconvolves time function TF out of records
%     in SEIZMO struct DATA using spectral division.  TF must be a numeric
%     vector or a cell array with one numeric vector per record in DATA.
%     The sample spacing of the time function is assumed to be the same as
%     the record it is deconvolved from.  Records in DATA must be evenly
%     spaced and Timeseries or XY datatype.  Note that if TF has spectral
%     amplitudes == 0, then the deconvolution becomes unstable (see the
%     next usage description for a more stable operation).
%
%     DATA=DECONVOLVE(DATA,TF,DELAY) specifies the delay of the time
%     function TF in samples.  By default DELAY is 0, which is the simple
%     causal convolution case.  Setting DELAY to a negative number
%     indicates the convolution was an acausal operation.
%
%     DATA=DECONVOLVE(DATA,TF,DELAY,H2O) adds factor H2O to the spectral
%     amplitudes of TF to avoid division by zero.  This is commonly known
%     as setting the waterlevel.  The best value for H2O varies and often
%     requires some inspection.  The default value is 0 and is the unstable
%     case.  H2O must be a positive scalar or an array of positive values,
%     with 1 value per record in DATA.  Typical values are 1e-4 to 1e-1.
%     Smaller values have less impact on the result while larger values
%     sacrifice accuracy for stability.  If the records in DATA require too
%     high of an H2O value, try the next option to exclude the unstable
%     frequencies.
%
%     DATA=DECONVOLVE(DATA,TF,DELAY,H2O,FRANGE) specifies the frequency
%     range for the spectral division.  Frequencies outside this range are
%     set to 0.  FRANGE must be 1x2 or Nx2 where N is the number of records
%     in DATA.  By default FRANGE is [0 NYQ] which includes all frequencies
%     (NYQ is the nyquist frequency for each record).  Use [] to get the
%     default range.
%
%     DATA=DECONVOLVE(DATA,TF,DELAY,H2O,FRANGE,ZI) removes convolution
%     initial conditions ZI from records in DATA before deconvolution.
%     This will remove the energy gained from adjacent data that can not be
%     accounted for otherwise.  See CONVOLVE for details.
%
%     DATA=DECONVOLVE(DATA,TF,DELAY,H2O,FRANGE,ZI,ZF) attaches convolution
%     final conditions ZF to the records in DATA.  This will account for
%     energy lost through the convolution operation to points extending
%     outside a record.  See CONVOLVE for details.
%
%    Notes:
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Convolve and then deconvolve a dataset with a waterlevel of 0.0001:
%     [data1,zf]=convolve(data,gausswin(13),-6);
%     data2=deconvolve(data1,gausswin(13),-6,0.0001,[],[],zf);
%
%    See also: CONVOLVE, IIRFILTER, CORRELATE

%     Version History:
%        Oct. 12, 2009 - initial version
%        Oct. 17, 2009 - added frange option
%        Oct. 21, 2009 - forced dim arg in fft/ifft
%        Oct. 28, 2009 - worked out (a)causal stuff, adds delay option,
%                        doubles initial conditions and final conditions
%        Jan. 27, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        fixed offset bug
%        Feb.  2, 2010 - versioninfo caching
%        Aug. 19, 2010 - removed ifft symmetric flag, real conversion
%        Feb. 11, 2011 - dropped versioninfo caching
%        Jan. 28, 2012 - doc update, better checkheader usage
%        May  30, 2012 - pow2pad=0 by default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt deconvolution
try
    % check headers
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % get header info
    [npts,ncmp,delta]=getheader(data,'npts','ncmp','delta');
    nyq=1./(2*delta);

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
        end
        ntf=numel(tf);
        tf={double(tf(:))};
        tfidx=ones(nrecs,1);
    end

    % check delay
    if(nargin<3 || isempty(delay)); delay=0; end
    if(~isreal(delay) || ~any(numel(delay)==[1 nrecs]) ...
            || ~all(delay==fix(delay)))
        error('seizmo:deconvolve:badDelay',...
            ['DELAY must be a scalar integer or an array\n' ...
            'of integers, one per record in DATA!']);
    end
    if(isscalar(delay)); delay=delay(ones(nrecs,1),1); end
    offset=max(0,delay);

    % check water level
    if(nargin<4 || isempty(h2o))
        h2o=0;
    elseif(~isreal(h2o) || (~isscalar(h2o) && numel(h2o)~=nrecs))
        error('seizmo:deconvolve:badH2O',...
            'H2O must be a scalar real or an array w/ 1 element/record!');
    elseif(any(h2o<0))
        error('seizmo:deconvolve:badH2O','H2O must be positive!');
    end
    if(isscalar(h2o)); h2o=h2o(ones(nrecs,1),1); end

    % check frequency range
    if(nargin<5 || isempty(frange))
        frange=[zeros(nrecs,1) nyq];
    elseif(~isnumeric(frange) || (~isequal(size(frange),[1 2]) ...
            && ~isequal(size(frange),[nrecs 2])))
        error('seizmo:deconvolve:badFRANGE',...
            'FRANGE must be a numeric array of size 1x2 or Nx2!')
    end
    if(size(frange,1)==1); frange=frange(ones(nrecs,1),:); end

    % get expected lengths
    nfz=max(0,ntf(tfidx)+delay-1);
    nbz=abs(min(0,delay));

    % check initial conditions
    haszi=nargin>5 && ~isempty(zi);
    if(haszi)
        % expand single columns
        if(iscell(zi) && size(zi,2)==1); zi{end,2}=[]; end
        if(~iscell(zi) || ~isequal(size(zi),[nrecs 2]))
            error('seizmo:deconvolve:badZI',...
                ['ZI must be a Nx2 cell array where\n'...
                'N is the number of records in DATA!']);
        end
        % allow empty if empty desired
        if(~all(cellfun('isreal',zi(:))) ...
                || ~all(...
                (cellfun('prodofsize',zi(:,1))==0 & nfz.*ncmp==0) ...
                | sum(([cellfun('size',zi(:,1),1) ...
                cellfun('size',zi(:,1),2)]==[nfz ncmp]),2)==2) ...
                || ~all(...
                (cellfun('prodofsize',zi(:,2))==0 & nbz.*ncmp==0) ...
                | sum(([cellfun('size',zi(:,2),1) ...
                cellfun('size',zi(:,2),2)]==[nbz ncmp]),2))==2)
            error('seizmo:deconvolve:badZI',...
                'ZI elements are improperly sized!');
        end
        for i=1:nrecs
            zi{i,1}=double(zi{i,1});
            zi{i,2}=double(zi{i,2});
        end
    end

    % check final conditions
    haszf=nargin>6 && ~isempty(zf);
    if(haszf)
        % expand single columns
        if(iscell(zf) && size(zf,2)==1); zf{end,2}=[]; end
        if(~iscell(zf) || ~isequal(size(zf),[nrecs 2]))
            error('seizmo:deconvolve:badZF',...
                ['ZF must be a Nx2 cell array where\n'...
                'N is the number of records in DATA!']);
        end
        % allow empty if empty desired
        if(~all(cellfun('isreal',zf(:))) ...
                || ~all(...
                (cellfun('prodofsize',zf(:,1))==0 & nfz.*ncmp==0) ...
                | sum(([cellfun('size',zf(:,1),1) ...
                cellfun('size',zf(:,1),2)]==[nfz ncmp]),2)==2) ...
                || ~all(...
                (cellfun('prodofsize',zf(:,2))==0 & nbz.*ncmp==0) ...
                | sum(([cellfun('size',zf(:,2),1) ...
                cellfun('size',zf(:,2),2)]==[nbz ncmp]),2))==2)
            error('seizmo:deconvolve:badZF',...
                'ZF elements are improperly sized!');
        end
        for i=1:nrecs
            zf{i,1}=double(zf{i,1});
            zf{i,2}=double(zf{i,2});
        end
    end
    
    % detail message
    if(verbose)
        disp('Deconvolving Time Function(s) from Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop through records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue; 
        end

        % get data class, convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % approach differs slightly for case with final conditions
        if(haszf)
            % attach final conditions
            data(i).dep=[zf{i,2}; data(i).dep; zf{i,1}];
        else
            % pad w/ zeros
            data(i).dep=[zeros(nbz(i),ncmp(i)); ...
                data(i).dep; zeros(nfz(i),ncmp(i))];
        end

        % start/stop indices of record
        start=nbz(i)+1;
        stop=nbz(i)+npts(i);

        % now remove initial conditions
        if(haszi)
            data(i).dep(start:(start+nfz(i)-1),:)=...
                data(i).dep(start:(start+nfz(i)-1),:)-zi{i,1};
            data(i).dep((stop-nbz(i)+1):stop,:)=...
                data(i).dep((stop-nbz(i)+1):stop,:)-zi{i,2};
        end

        % waterleveled spectral division
        nspts=2^nextpow2(npts(i)+ntf(tfidx(i))-1);
        sdelta=2*nyq(i)./nspts;
        freq=abs([0:(nspts/2) (1-nspts/2):-1].*sdelta);
        good=freq>=frange(i,1) & freq<=frange(i,2);
        tmp1=fft(data(i).dep,nspts,1);
        tmp2=fft(tf{tfidx(i)}(:,ones(1,ncmp(i))),nspts,1);
        tmp2=(abs(tmp2)+h2o(i)).*exp(1j*angle(tmp2));
        tmp=complex(zeros(nspts,ncmp(i)));
        tmp(good,:)=tmp1(good,:)./tmp2(good,:);
        tmp=real(ifft(tmp,[],1));
        data(i).dep=tmp(offset(i)+(1:npts(i)),:);

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));

        % restore class
        data(i).dep=oclass(data(i).dep);
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,...
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
