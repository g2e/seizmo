function [data]=taper(data,type,width,option)
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
%   Usage:  [data]=taper(data,type,width,option)
%
%   Examples:
%    Taper data with a gaussian that is applied to the first and last 10th
%    of the record with an option selected that makes the gaussian taper
%    represent a gaussian curve from peak out to 4 standard deviations
%      data=taper(data,'gausswin',0.1,4);
%
%    Default tapering - blackman-harris taper applied to the first and last
%    20th of the record
%      data=taper(data);
%
%    See also: rmean, rtrend, rdrift

% check input
error(nargchk(1,4,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% defaults
if(nargin<2 || isempty(type)); type='blackmanharris'; end
if(nargin<3 || isempty(width)); width=0.05; end

% check width
if(width>1); error('taper halfwidth far too big - use 0 to 0.5'); end

% make function handle
type=str2func(type);

% header info
[leven,iftype]=gh(data,'leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% work through each file
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % check if unevenly spaced
    if(leven~=h(v).true)
        warning('SAClab:illegalOperation','Illegal operation on unevenly spaced record %d',i);
        continue;
    end
    
    % check for unsupported filetypes
    if(iftype(i)==h(v).enum(1).val.ixyz)
        warning('SAClab:illegalFiletype','Illegal operation on xyz file');
        continue;
    end
    
    % record size, taper size
    [npts,ncmp]=size(data(i).x);
    nwidth=ceil(width*npts);
    
    % make taper
    if(nargin==4 && ~isempty(option)); taperedge=window(type,2*nwidth,option);
    else taperedge=window(type,2*nwidth); end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % apply taper halfwidths separately
    data(i).x(1:nwidth,:)=...
        data(i).x(1:nwidth,:).*taperedge(1:nwidth,ones(ncmp,1));
    data(i).x(end-nwidth+1:end,:)=...
        data(i).x(end-nwidth+1:end,:).*taperedge(end-nwidth+1:end,ones(ncmp,1));
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end