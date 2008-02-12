function [data]=rslope(data)
%RTREND    Romove slope from SAClab data records
%
%    Description: Removes the slope from records while preserving the data
%     mean (as best as possible).  Uses Matlab functions detrend or polyfit
%     and polyval.
%
%    Usage:  [data]=rslope(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: rmean, rdrift

% check input
error(nargchk(1,1,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% grab timing info
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end
leven=gh(data,'leven');

% remove trend and update header
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % act based on data spacing
    if(leven(i)==h(v).true)
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=detrend(data(i).x(:,j))+mean(data(i).x(:,j));
        end
    else
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=data(i).x(:,j) ...
                -polyval(polyfit(double(data(i).t),data(i).x(:,j),1),...
                double(data(i).t))+mean(data(i).x(:,j));
        end
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end
