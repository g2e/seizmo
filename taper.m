function [data]=taper(data,width,type,option)
%TAPER   Taper SEIZMO records
%
%    Usage:    data=taper(data)
%              data=taper(data,width)
%              data=taper(data,width,type)
%              data=taper(data,width,type,option)
%
%    Description: TAPER(DATA) tapers data records with a Hann taper set to
%     vary from 0 to 1 over 5% of every records' length on each end.  This
%     matches SAC's default taper command.
%
%     TAPER(DATA,WIDTH) allows the taper width to be set.  WIDTH should be
%     a number anywhere from 0.0 (no taper) to 0.5 (taper entire record).
%     Two numbers may be given to apply a different width taper to each
%     end (first number gives leading taper).
%
%     TAPER(DATA,WIDTH,TYPE) allows the taper type to be changed.  TYPE is
%     a string that must be one of the following:
%
%       TYPE string     Formal window/taper name       Option
%       %%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%       %%%%%%%%%%%%%%%%%%%%
%       barthannwin     Modified Bartlett-Hann         
%       bartlett        Bartlett                       
%       blackman        Blackman                       'periodic|symmetric'
%       blackmanharris  Minimum 4-term Blackman-Harris 
%       bohmanwin       Bohman                         
%       chebwin         Chebyshev                       sidelobe_atten*
%       flattopwin      Flat Top weighted              'periodic|symmetric'
%       gausswin        Gaussian                        num_std_dev
%       hamming         Hamming                        'periodic|symmetric'
%       hann            Hann (Hanning)                 'periodic|symmetric'
%       kaiser          Kaiser                          beta_parameter*
%       nuttallwin      Nuttall-defined Blackman-Harris
%       parzenwin       Parzen (de la Valle-Poussin)
%       rectwin         Rectangular (no taper)
%       triang          Triangular 
%       tukeywin        Tukey (tapered cosine)          taper_ratio*
%
%     More information on each taper can be found with 'help <type_string>'
%     and 'doc <type_string>' where <type_string> should be replaced with
%     the taper's above TYPE string.
%
%     TAPER(DATA,WIDTH,TYPE,OPTION) sets a taper parameter to OPTION.  For
%     taper types 'chebwin', 'kaiser' and 'tukeywin' the taper parameter is
%     required.  Use 'help <type_string>' and 'doc <type_string>', where 
%     <type_string> should be replaced with the taper's above TYPE string,
%     for specifics on each taper's parameter.
%
%    Notes:
%     - uses Matlab's Signal Processing Toolbox WINDOW function
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Taper data with a gaussian that is applied to the first and last 10th
%     of the record with the taper forced to represent a gaussian curve 
%     from peak out to 4 standard deviations:
%      data=taper(data,0.1,'gausswin',4);
%
%    See also: removemean, removetrend

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2009 at 20:00 GMT

% check input
msg=nargchk(1,4,nargin);
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
if(nargin<3 || isempty(type)); type='hann'; end
if(nargin<2 || isempty(width)); width=[0.05 0.05]; end

% check width
if(numel(width)==1)
    width=[width width];
elseif(numel(width)~=2)
    error('seizmo:taper:badInput','Too many taper halfwidth parameters!');
end
if(any(width<0))
    error('seizmo:taper:badInput','Taper halfwidths must be >=0!');
end

% make function handle
type=str2func(type);

% header info
leven=getlgc(data,'leven');
iftype=getenumdesc(data,'iftype');

% check for unsupported filetypes
if(any(strcmpi(iftype,'General XYZ (3-D) file')))
    error('seizmo:taper:illegalFiletype',...
        'Illegal operation on xyz file');
end

% number of records
nrecs=numel(data);

% work through each file
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % record size, taper size
    [npts,ncmp]=size(data(i).dep);
    nwidth=ceil(width*npts);
    
    % unevenly spaced
    if(strcmp(leven(i),'false'))
        % make tapers
        if(nargin==4 && ~isempty(option))
            taperedge1=window(type,2*nwidth(1),option);
            taperedge2=window(type,2*nwidth(2),option);
        else
            taperedge1=window(type,2*nwidth(1));
            taperedge2=window(type,2*nwidth(2));
        end
        
        % interpolate
        last1=find(data(i).ind<(data(i).ind(1)+...
            (data(i).ind(end)-data(i).ind(1))*width(1)),1,'last');
        last2=find(data(i).ind>(data(i).ind(end)-...
            (data(i).ind(end)-data(i).ind(1))*width(2)),1,'first');
        even_times=linspace(data(i).ind(1),data(i).ind(end),npts);
        taper1=interp1(even_times(1:nwidth(1)),taperedge1(1:nwidth(1)),...
            data(i).ind(1:last1),'spline');
        taper2=interp1(even_times((end-nwidth(2)+1):end),...
            taperedge2((nwidth(2)+1):end),data(i).ind(last2:end),'spline');
        
        % apply taper halfwidths separately
        data(i).dep(1:last1,:)=...
            data(i).dep(1:last1,:).*taper1(:,ones(ncmp,1));
        data(i).dep(last2:end,:)=...
            data(i).dep(last2:end,:).*taper2(:,ones(ncmp,1));
    % evenly spaced
    else
        % time series and general xy records
        if(strcmp(iftype(i),'itime') || strcmp(iftype(i),'ixy'))
            % make tapers
            if(nargin==4 && ~isempty(option))
                taperedge1=window(type,2*nwidth(1),option);
                taperedge2=window(type,2*nwidth(2),option);
            else
                taperedge1=window(type,2*nwidth(1));
                taperedge2=window(type,2*nwidth(2));
            end
            
            % adjust for npts
            nwidth=min(npts,nwidth);
            
            % apply taper halfwidths separately
            data(i).dep(1:nwidth(1),:)=data(i).dep(1:nwidth(1),:)...
                .*taperedge1(1:nwidth(1),ones(ncmp,1));
            data(i).dep((end-nwidth(2)+1):end,:)=...
                data(i).dep((end-nwidth(2)+1):end,:)...
                .*taperedge2((end-nwidth(2)+1):end,ones(ncmp,1));
        else % spectral
            % fix nwidth
            nwidth=ceil(width*npts/2);
            
            % make tapers
            if(nargin==4 && ~isempty(option))
                taperedge1=window(type,2*nwidth(1),option);
                taperedge2=window(type,2*nwidth(2),option);
            else
                taperedge1=window(type,2*nwidth(1));
                taperedge2=window(type,2*nwidth(2));
            end
            
            % adjust for npts
            nwidth=min(npts/2,nwidth);
            
            % apply taper halfwidths separately
            data(i).dep(1:nwidth(1),:)=data(i).dep(1:nwidth(1),:)...
                .*taperedge1(1:nwidth(1),ones(ncmp,1));
            data(i).dep((end-nwidth(1)+2):end,:)=...
                data(i).dep((end-nwidth(1)+2):end,:)...
                .*taperedge1((end-nwidth(1)+2):end,ones(ncmp,1));
            nyqpt=npts/2+1;
            data(i).dep((nyqpt-nwidth(2)+1):nyqpt,:)=...
                data(i).dep((nyqpt-nwidth(2)+1):nyqpt,:)...
                .*taperedge2((end-nwidth(2)+1):end,ones(ncmp,1));
            data(i).dep((nyqpt+1):(nyqpt+nwidth(2)-1),:)=...
                data(i).dep((nyqpt+1):(nyqpt+nwidth(2)-1),:)...
                .*taperedge2(2:nwidth(2),ones(ncmp,1));
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
