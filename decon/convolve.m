function [data,zf]=convolve(data,tf,delay,zi)
%CONVOLVE    Convolve SEIZMO records with a time function
%
%    Usage:    data=convolve(data,tf)
%              [data,zf]=convolve(data,tf)
%              [data,zf]=convolve(data,tf,delay)
%              [data,zf]=convolve(data,tf,delay,zi)
%
%    Description:
%     DATA=CONVOLVE(DATA,TF) convolves a time function TF onto records in
%     SEIZMO struct DATA.  TF must be a numeric vector or a cell array with
%     as many elements as records in DATA with each element being a numeric
%     vector (allows for each record to have a different time function
%     applied).  The sample spacing of the time function is assumed to be
%     the same as that of each record.  Convolution is done in double
%     precision in the time domain using FILTER.  The records are then
%     converted back to their original class.
%
%     [DATA,ZF]=CONVOLVE(DATA,TF) returns the final conditions of the
%     convolution in ZF.  ZF is a Nx2 cell array with each row
%     corresponding to an individual record in DATA.  Column 1 of ZF gives
%     the forward (causal) final conditions while column 2 contains the
%     backward (acausal) final conditions.  So ZF(3,:) corresponds
%     to DATA(3).  Final conditions are particularly useful to communicate
%     convolution info between sequential segments of a seismogram by
%     passing the info as initial conditions to adjacent segments.
%
%     [DATA,ZF]=CONVOLVE(DATA,TF,DELAY) specifies the delay of the time
%     function TF in samples.  By default DELAY is 0, which is the simple
%     causal convolution case.  Setting DELAY to a negative number
%     indicates the convolution is an acausal operation.  See the Examples
%     section below for an example of using the DELAY argument.
%
%     [DATA,ZF]=CONVOLVE(DATA,TF,DELAY,ZI) applies initial conditions ZI to
%     the convolution.  This is useful when you have sequential records
%     (see the description of output ZF for more details).
%
%    Notes:
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Convolve a record with a centered 11-point gaussian:
%     [data1,zf]=convolve(data(1),gausswin(11),-5);
%
%     % And plot an overlay to confirm (attaching final conditions too):
%     plot2([data(1); ...
%         attach(attach(data1,'ending',zf{1,1}),'beginning',zf{1,2})])
%
%    See also: DECONVOLVE, ATTACH, IIRFILTER, CORRELATE

%     Version History:
%        Oct.  8, 2009 - initial version
%        Oct. 10, 2009 - LEVEN and IFTYPE checks
%        Oct. 28, 2009 - worked out (a)causal stuff, adds delay option,
%                        doubles initial conditions and final conditions,
%                        fixed a couple bugs (dataless, zi applied wrong)
%        Jan. 27, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  2, 2010 - versioninfo caching
%        Mar.  8, 2010 - dropped versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix
%        Jan. 28, 2012 - doc update, better checkheader usage
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt convolution
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
    ncmp=getheader(data,'ncmp');

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

    % check delay
    if(nargin<3 || isempty(delay)); delay=0; end
    if(~isreal(delay) || ~any(numel(delay)==[1 nrecs]) ...
            || ~all(delay==fix(delay)))
        error('seizmo:convolve:badDelay',...
            ['DELAY must be a scalar integer or an array\n' ...
            'of integers, one per record in DATA!']);
    end
    if(isscalar(delay)); delay=delay(ones(nrecs,1),1); end

    % get expected lengths
    nfz=max(0,ntf(tfidx)+delay-1);
    nbz=abs(min(0,delay));

    % check initial conditions
    haszi=nargin>3 && ~isempty(zi);
    if(haszi)
        % expand single columns
        if(iscell(zi) && size(zi,2)==1); zi{end,2}=[]; end
        if(~iscell(zi) || ~isequal(size(zi),[nrecs 2]))
            error('seizmo:convolve:badZI',...
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
            error('seizmo:convolve:badZI',...
                'ZI elements are improperly sized!');
        end
        for i=1:nrecs
            zi{i,1}=double(zi{i,1});
            zi{i,2}=double(zi{i,2});
        end
    end
    
    % detail message
    if(verbose)
        disp('Convolving Time Function(s) on Record(s)');
        print_time_left(0,nrecs);
    end

    % loop through records
    zf=cell(nrecs,2);
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % simple ascii illustration of how this looks when
        % the record is shorter than the neg/pos delays of
        % the time function (zi alter the zf):
        %
        % dep in:               ----
        %
        % dep out:              ----
        % zf  for:                  -------
        % zf  bak:     ---------
        % zi  for:              -------
        % zi  bak:         ---------

        % get data class, convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % skip dataless (do convert class back)
        if(isempty(data(i).dep))
            % note zf==zi in this case
            data(i).dep=oclass(data(i).dep);
            zf{i,1}=oclass(zi{i,1});
            zf{i,2}=oclass(zi{i,2});
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % convolve (assumes simple causal case)
        [data(i).dep,zf{i,1}]=filter(tf{tfidx(i)},1,data(i).dep,[],1);

        % combined convolved record with final conditions (padding too)
        data(i).dep=...
            [zeros(max(0,delay(i)),ncmp(i)); data(i).dep; zf{i,1}; ...
            zeros(abs(min(ntf(tfidx(i))-1-nbz(i),0)),ncmp(i))];

        % start/stop indice of record
        start=nbz(i)+1;
        stop=size(data(i).dep,1)-nfz(i);

        % now apply initial conditions ontop of everything
        if(haszi)
            data(i).dep(start:(start+nfz(i)-1),:)=...
                data(i).dep(start:(start+nfz(i)-1),:)+zi{i,1};
            data(i).dep((stop-nbz(i)+1):stop,:)=...
                data(i).dep((stop-nbz(i)+1):stop,:)+zi{i,2};
        end

        % extract final conditions (leaving just the record)
        zf{i,1}=data(i).dep((stop+1):end,:);
        zf{i,2}=data(i).dep(1:nbz(i),:);
        data(i).dep([1:nbz(i) (stop+1):end],:)=[];

        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));

        % restore class
        data(i).dep=oclass(data(i).dep);
        zf{i,1}=oclass(zf{i,1});
        zf{i,2}=oclass(zf{i,2});
        
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
