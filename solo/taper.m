function [data]=taper(data,w,o,type,opt)
%TAPER   Taper SEIZMO records
%
%    Usage:    data=taper(data)
%              data=taper(data,width)
%              data=taper(data,width,offset)
%              data=taper(data,width,offset,type)
%              data=taper(data,width,offset,type,option)
%
%    Description:
%     TAPER(DATA) tapers data records with a Hann taper set to vary from 0
%     to 1 over 5% of every records' length on each end.  This matches
%     SAC's default taper command.
%
%     TAPER(DATA,WIDTH) allows the taper width to be set.  WIDTH should be
%     a number anywhere from 0.0 (no taper) to 0.5 (taper from the edge to
%     the center of the record).  The tapering is applied symmetrically (as
%     in the taper is applied to both ends of the records).  To specify
%     different leading and trailing tapers set WIDTH to [WIDTH1 WIDTH2],
%     where WIDTH1 gives the width of the leading taper and WIDTH2
%     indicates the width of the trailing taper.  The default width is
%     0.05.
%
%     TAPER(DATA,WIDTH,OFFSET) specifies the offset of the tapers from the
%     start/end of the records as a fraction of the records width.  Any
%     points in the record nearer to the start/end of the record than
%     OFFSET are set to 0.  So an offset of 0.1 for a taper of width 0.2
%     will zero out the first and last tenth of the record and taper from
%     0.1 to 0.3 and from 0.7 to 0.9.  The default offset is 0.
%
%     TAPER(DATA,WIDTH,OFFSET,TYPE) allows the taper type to be changed.
%     TYPE is a string that must be one of the following:
%
%       TYPE string     Formal window/taper name       Option
%       %%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%       %%%%%%%%%%%%%%%%%%%%
%       barthannwin     Modified Bartlett-Hann         
%       bartlett        Bartlett                       
%       blackman        Blackman                       
%       blackmanharris  Minimum 4-term Blackman-Harris 
%       bohmanwin       Bohman                         
%       chebwin         Chebyshev                      sidelobe_atten (100)
%       flattopwin      Flat Top weighted              
%       gausswin        Gaussian                       num_std_dev (2.5)
%       hamming         Hamming                        
%       hann            Hann (Hanning)                 
%       kaiser          Kaiser                         beta_parameter (0.5)
%       nuttallwin      Nuttall-defined Blackman-Harris
%       parzenwin       Parzen (de la Valle-Poussin)   
%       rectwin         Rectangular (no taper)         
%       triang          Triangular                     
%       tukeywin        Tukey (tapered cosine)         taper_ratio (0.5)
%
%     More information on each taper can be found with 'help <type_string>'
%     and 'doc <type_string>' where <type_string> should be replaced with
%     the taper's above TYPE string.  The default type is 'hann'.
%
%     TAPER(DATA,WIDTH,OFFSET,TYPE,OPTION) sets a taper's parameter to
%     OPTION.  Use 'help <type_string>' and 'doc <type_string>', where 
%     <type_string> should be replaced with the taper's above TYPE string,
%     for specifics on each taper's parameter.  Setting OPTION to [] or nan
%     uses the default value.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Taper data with a gaussian that is applied to the first and last
%     % 10th of the record with the taper forced to represent a gaussian
%     % curve from peak out to 4 standard deviations:
%     data=taper(data,0.1,0,'gausswin',4);
%
%     % Now do a similar taper but start the leading taper after the first
%     % third and the trailing taper before the final eighth:
%     data=taper(data,0.1,[1/3 0.125],'gausswin',4)
%
%    See also: REMOVEMEAN, REMOVETREND, REMOVEPOLYNOMIAL

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Feb. 12, 2008 - rename from SACTAPER to TAPER, handle string input
%        Feb. 16, 2008 - class support, multi-component support
%        Feb. 25, 2008 - support for unevenly sampled records, better
%                        checks, input arguments order changed
%        Mar.  4, 2008 - minor doc update, major code cleaning
%        Apr. 17, 2008 - minor doc update
%        May  12, 2008 - dep* formula fix
%        June 12, 2008 - doc update, made 'hann' default to match SAC
%        Nov. 22, 2008 - doc update, history fix, update for new name
%                        schema, handle widths>1, one changeheader call,
%                        error on xyz file, better checking
%        Dec. 12, 2008 - doc update
%        Apr. 22, 2009 - changed interpolation method for building tapers
%                        related to unevenly sampled records
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 11, 2009 - special handling of spectral records (4 tapers)
%        June 29, 2009 - fix iftype bug
%        Sep. 23, 2009 - major code revision, changed input args, allow for
%                        1 entry per record, more checks, use new taperfun
%                        function to handle taper building
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%        Jan. 30, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  2, 2010 - versioninfo caching
%        Aug. 25, 2010 - drop versioninfo caching, nargchk fix
%        Mar. 13, 2012 - doc update, seizmocheck fix, better checkheader
%                        usage, use getheader improvements
%        Feb. 14, 2013 - using strcmpi for consistency
%        Feb. 26, 2013 - bugfix: workaround precision issues for 0 width
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 21:00 GMT

% todo:

% check input
error(nargchk(1,5,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt taper
try
    % check headers
    data=checkheader(data,...
        'XYZ_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % defaults
    if(nargin<5 || isempty(opt)); opt=nan; end
    if(nargin<4 || isempty(type)); type='hann'; end
    if(nargin<3 || isempty(o)); o=[0 0]; end
    if(nargin<2 || isempty(w)); w=[0.05 0.05]; end

    % check width (class, # of rows, # of columns)
    if(~isreal(w))
        error('seizmo:taper:badInput','WIDTH must be real!');
    end
    sz=size(w);
    if(sz(1)==1)
        w=w(ones(nrecs,1),:);
    elseif(sz(1)~=nrecs)
        error('seizmo:taper:badInput',...
            'WIDTH must have 1 row or 1 row per record!');
    end
    if(sz(2)==1)
        w=[w w];
    elseif(sz(2)~=2)
        error('seizmo:taper:badInput',...
            'WIDTH must have 1 or 2 columns!');
    end

    % fix negative width
    if(any(w<0)); w(w<0)=0; end

    % check offsets (class, # of rows, # of columns)
    if(~isreal(o))
        error('seizmo:taper:badInput','OFFSET must be real!');
    end
    sz=size(o);
    if(sz(1)==1)
        o=o(ones(nrecs,1),:);
    elseif(sz(1)~=nrecs)
        error('seizmo:taper:badInput',...
            'OFFSET must have 1 row or 1 row per record!');
    end
    if(sz(2)==1)
        o=[o o];
    elseif(sz(2)~=2)
        error('seizmo:taper:badInput',...
            'OFFSET must have 1 or 2 columns!');
    end

    % convert width to limit
    lim1=[o(:,1) o(:,1)+w(:,1)];
    lim2=[o(:,2) o(:,2)+w(:,2)];

    % check type (class, # of rows, # of columns)
    if(ischar(type)); type=cellstr(type); end
    sz=size(type);
    if(~iscellstr(type))
        error('seizmo:taper:badInput',...
            'TYPE must be a char or cellstr array!');
    end
    if(sz(1)==1)
        type=type(ones(nrecs,1),:);
    elseif(sz(1)~=nrecs)
        error('seizmo:taper:badInput',...
            'TYPE must have 1 row or 1 row per record!');
    end
    if(sz(2)==1)
        type=[type type];
    elseif(sz(2)~=2)
        error('seizmo:taper:badInput',...
            'TYPE cell array must have 1 or 2 columns!');
    end

    % check option (class, # of rows, # of columns)
    if(~isreal(opt))
        error('seizmo:taper:badInput','OPTION must be real!');
    end
    sz=size(opt);
    if(sz(1)==1)
        opt=opt(ones(nrecs,1),:);
    elseif(sz(1)~=nrecs)
        error('seizmo:taper:badInput',...
            'OPTION must have 1 row or 1 row per record!');
    end
    if(sz(2)==1)
        opt=[opt opt];
    elseif(sz(2)~=2)
        error('seizmo:taper:badInput',...
            'OPTION must have 1 or 2 columns!');
    end

    % header info
    [b,e,delta,npts,ncmp,leven,iftype]=getheader(data,...
        'b','e','delta','npts','ncmp','leven lgc','iftype id');
    
    % detail message
    if(verbose)
        disp('Tapering Record(s)');
        print_time_left(0,nrecs);
    end

    % work through each file
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % unevenly spaced
        if(strcmpi(leven(i),'false'))
            % get normalized time
            times=[0 (data(i).ind(2:end-1)'-b(i))/(e(i)-b(i)) 1];
            if(npts(i)<2); times=0; end

            % get tapers
            taper1=taperfun(type{i,1},times,lim1(i,:),opt(i,1));
            taper2=taperfun(type{i,2},1-times,lim2(i,:),opt(i,2));
            taper1=taper1(:); taper2=taper2(:);

            % apply tapers
            data(i).dep=data(i).dep.*taper1(:,ones(1,ncmp(i))) ...
                .*taper2(:,ones(1,ncmp(i)));
        % evenly spaced
        else
            % time series and general xy records
            if(strcmpi(iftype(i),'itime') || strcmpi(iftype(i),'ixy'))
                % get normalized time
                times=delta(i)/(e(i)-b(i));
                times=[0 times:times:times*(npts(i)-2) 1];
                if(npts(i)<2); times=0; end

                % get tapers
                taper1=taperfun(type{i,1},times,lim1(i,:),opt(i,1));
                taper2=taperfun(type{i,2},1-times,lim2(i,:),opt(i,2));
                taper1=taper1(:); taper2=taper2(:);

                % apply tapers
                data(i).dep=data(i).dep.*taper1(:,ones(1,ncmp(i))) ...
                    .*taper2(:,ones(1,ncmp(i)));
            else % spectral
                % get normalized frequency
                freq=delta(i)/e(i);
                freq=[0:freq:freq*(npts(i)/2-1) 1 ...
                    freq*(npts(i)/2-1):-freq:freq];
                if(npts(i)<2); freq=0; end

                % get tapers
                taper1=taperfun(type{i,1},freq,lim1(i,:),opt(i,1));
                taper2=taperfun(type{i,2},1-freq,lim2(i,:),opt(i,2));
                taper1=taper1(:); taper2=taper2(:);

                % apply tapers
                data(i).dep=data(i).dep.*taper1(:,ones(1,2*ncmp(i))) ...
                    .*taper2(:,ones(1,2*ncmp(i)));
            end
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % get dep* values
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
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
