function [data]=taper(data,w,o,type,opt)
%TAPER   Taper SEIZMO records
%
%    Usage:    data=taper(data)
%              data=taper(data,width)
%              data=taper(data,width,offset)
%              data=taper(data,width,offset,type)
%              data=taper(data,width,offset,type,option)
%
%    Description: TAPER(DATA) tapers data records with a Hann taper set to
%     vary from 0 to 1 over 5% of every records' length on each end.  This
%     matches SAC's default taper command.
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
%     Taper data with a gaussian that is applied to the first and last 10th
%     of the record with the taper forced to represent a gaussian curve 
%     from peak out to 4 standard deviations:
%      data=taper(data,0.1,0,'gausswin',4);
%
%     Now do a similar taper but start the leading taper after the first
%     third and the trailing taper before the final eighth:
%      data=taper(data,0.1,[1/3 0.125],'gausswin',4)
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 23, 2009 at 18:15 GMT

% todo:

% check input
msg=nargchk(1,5,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% defaults
if(nargin<5 || isempty(opt)); opt=nan; end
if(nargin<4 || isempty(type)); type='hann'; end
if(nargin<3 || isempty(o)); o=[0 0]; end
if(nargin<2 || isempty(w)); w=[0.05 0.05]; end

% number of records
nrecs=numel(data);

% check width (class, # of rows, # of columns)
if(~isreal(w)); error('seizmo:taper:badInput','WIDTH must be real!'); end
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
if(~isreal(o)); error('seizmo:taper:badInput','OFFSET must be real!'); end
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
        'OPTON must have 1 row or 1 row per record!');
end
if(sz(2)==1)
    opt=[opt opt];
elseif(sz(2)~=2)
    error('seizmo:taper:badInput',...
        'OPTION must have 1 or 2 columns!');
end

% header info
[b,e,delta,npts,ncmp]=getheader(data,'b','e','delta','npts','ncmp');
leven=getlgc(data,'leven');
iftype=getenumid(data,'iftype');

% check for unsupported filetypes
if(any(strcmpi(iftype,'ixyz')))
    error('seizmo:taper:illegalFiletype',...
        'Illegal operation on xyz file');
end

% work through each file
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % unevenly spaced
    if(strcmp(leven(i),'false'))
        % get normalized time
        times=(data(i).ind-b(i))/(e(i)-b(i));
        
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
        if(strcmp(iftype(i),'itime') || strcmp(iftype(i),'ixy'))
            % get normalized time
            times=(linspace(b(i),e(i),npts(i))-b(i))/(e(i)-b(i));
            
            % get tapers
            taper1=taperfun(type{i,1},times,lim1(i,:),opt(i,1));
            taper2=taperfun(type{i,2},1-times,lim2(i,:),opt(i,2));
            taper1=taper1(:); taper2=taper2(:);
            
            % apply tapers
            data(i).dep=data(i).dep.*taper1(:,ones(1,ncmp(i))) ...
                .*taper2(:,ones(1,ncmp(i)));
        else % spectral
            % get normalized frequency
            freq=[linspace(b(i),e(i),npts(i)/2+1) linspace( ...
                e(i)-delta(i),b(i)+delta(i),npts(i)/2-1)]/(e(i)-b(i));
            
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
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
