function [data]=taper(data,width,type,option)
%TAPER   Taper SAClab data records
%
%    Description: A taper is a monotonically varying function between zero
%     and one.  It is applied in a symmetric manner to the data such that 
%     the signal is zero for the first and last data points and increases 
%     smoothly to its original value at an interior point relative to each 
%     end.  This command allows SAClab to utilize the multitude of window
%     functions available in Matlab's Signal Processing Toolbox for tapers.
%
%     TAPER(DATA) tapers data records with a Hanning taper set to vary from
%     0 to 1 over 0.05 of every records' width on each end.  This matches
%     SAC's default taper command.
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
%     required.  Use Matlab's help for specifics on each taper's parameter.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%   
%    Usage: data=taper(data)
%           data=taper(data,width)
%           data=taper(data,width,type)
%           data=taper(data,width,type,option)
%
%    Examples:
%     Taper data with a gaussian that is applied to the first and last 10th
%     of the record with the taper forced to represent a gaussian curve 
%     from peak out to 4 standard deviations:
%      data=taper(data,0.1,'gausswin',4);
%
%    See also: rmean, rtrend

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Cleaned up documentation, made 'hann' default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 04:35 GMT

% check input
error(nargchk(1,4,nargin))

% check data structure
error(seischk(data,'x'))

% defaults
if(nargin<3 || isempty(type)); type='hann'; end
if(nargin<2 || isempty(width)); width=[0.05 0.05]; end

% check width
if(length(width)==1); width=[width width];
elseif(length(width)~=2)
    error('SAClab:taper:badInput','too many taper halfwidth parameters')
end
if(any(width>1))
    error('SAClab:taper:badInput',...
        'taper halfwidth far too big - use 0.0 to 0.5')
end

% make function handle
type=str2func(type);

% header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');

% check sample spacing logical
tru=strcmp(leven,'true');
fals=strcmp(leven,'false');
if(~all(tru | fals))
    error('SAClab:taper:levenBad',...
        'logical field leven needs to be set'); 
end

% work through each file
for i=1:length(data)
    % check for unsupported filetypes
    if(strcmp(iftype(i),'General XYZ (3-D) file'))
        warning('SAClab:taper:illegalFiletype',...
            'Illegal operation on xyz file');
        continue;
    end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % record size, taper size
    [npts,ncmp]=size(data(i).x);
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
        last1=find(data(i).t>data(i).t(1)+(data(i).t(end)-data(i).t(1))*width(1),1)-1;
        last2=find(data(i).t<data(i).t(end)-(data(i).t(end)-data(i).t(1))*width(2),1,'last')+1;
        even_times=linspace(data(i).t(1),data(i).t(end),npts);
        taper1=interp1(even_times(1:nwidth(1)),taperedge1(1:nwidth(1)),data(i).t(1:last1),'pchip');
        taper2=interp1(even_times(end-nwidth(2)+1:end),taperedge2(nwidth(2)+1:end),data(i).t(last2:end),'pchip');
        
        % apply taper halfwidths separately
        data(i).x(1:last1,:)=data(i).x(1:last1,:).*taper1(:,ones(ncmp,1));
        data(i).x(last2:end,:)=data(i).x(last2:end,:).*taper2(:,ones(ncmp,1));
    % evenly spaced
    else
        % make tapers
        if(nargin==4 && ~isempty(option))
            taperedge1=window(type,2*nwidth(1),option);
            taperedge2=window(type,2*nwidth(2),option);
        else
            taperedge1=window(type,2*nwidth(1));
            taperedge2=window(type,2*nwidth(2));
        end
        
        % apply taper halfwidths separately
        data(i).x(1:nwidth(1),:)=data(i).x(1:nwidth(1),:).*taperedge1(1:nwidth(1),ones(ncmp,1));
        data(i).x(end-nwidth(2)+1:end,:)=data(i).x(end-nwidth(2)+1:end,:).*taperedge2(nwidth(2)+1:end,ones(ncmp,1));
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
