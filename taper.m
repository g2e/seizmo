function [data]=taper(data,width,type,option)
%TAPER   Taper SAClab data records
%
%    Description: Tapers data records with the specified taper type and
%     width (fraction of total signal length).  The taper is applied 
%     symmetrically so that when the total width is less than 1, the taper 
%     leaves the middle portion untapered.  Optional inputs are the taper 
%     function type input as a string ('blackmanharris', 'hamming', 'hann', 
%     'gausswin', etc - look under Matlab's window function for all the 
%     types - default is 'blackmanharris'), the width parameter (ratio of 
%     taper halfwidth to the whole record length - 0 to 0.5 - default is 
%     0.05 which tapers the first and last 20th of the record), and an 
%     option to pass to Matlab's window function (for more precise control 
%     over certain taper types - read Matlab's window command for info).
%
%   Usage:  [data]=taper(data,width,type,option)
%
%   Examples:
%    Taper data with a gaussian that is applied to the first and last 10th
%    of the record with an option selected that makes the gaussian taper
%    represent a gaussian curve from peak out to 4 standard deviations
%      data=taper(data,0.1,'gausswin',4);
%
%    Default tapering - blackman-harris taper applied to the first and last
%    20th of the record
%      data=taper(data);
%
%    See also: rmean, rtrend, rdrift

% check input
error(nargchk(1,4,nargin))

% check data structure
if(~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% defaults
if(nargin<3 || isempty(type)); type='blackmanharris'; end
if(nargin<2 || isempty(width)); width=[0.05 0.05]; end

% check width
if(length(width)==1); width=[width width];
elseif(length(width)~=2); error('too many taper halfwidth parameters'); end
if(any(width>1)); error('taper halfwidth far too big - use 0 to 0.5'); end

% make function handle
type=str2func(type);

% header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');

% work through each file
for i=1:length(data)
    % check for unsupported filetypes
    if(strcmp(iftype(i),'General XYZ (3-D) file'))
        warning('SAClab:illegalFiletype','Illegal operation on xyz file');
        continue;
    end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % record size, taper size
    [npts,ncmp]=size(data(i).x);
    nwidth=ceil(width*npts);
    
    % check if unevenly spaced
    if(strcmp(leven(i),'false'))
        % make tapers
        if(nargin==4 && ~isempty(option)); taperedge1=window(type,2*nwidth(1),option);
        else taperedge1=window(type,2*nwidth(1)); end
        if(nargin==4 && ~isempty(option)); taperedge2=window(type,2*nwidth(2),option);
        else taperedge2=window(type,2*nwidth(2)); end
        
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
    elseif(strcmp(leven(i),'true'))
        % make taper
        if(nargin==4 && ~isempty(option)); taperedge1=window(type,2*nwidth(1),option);
        else taperedge1=window(type,2*nwidth(1)); end
        if(nargin==4 && ~isempty(option)); taperedge2=window(type,2*nwidth(2),option);
        else taperedge2=window(type,2*nwidth(2)); end
        
        % apply taper halfwidths separately
        data(i).x(1:nwidth(1),:)=data(i).x(1:nwidth(1),:).*taperedge1(1:nwidth(1),ones(ncmp,1));
        data(i).x(end-nwidth(2)+1:end,:)=data(i).x(end-nwidth(2)+1:end,:).*taperedge2(nwidth(2)+1:end,ones(ncmp,1));
    else
        error('sample spacing logical needs to be set for record %d',i)
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end