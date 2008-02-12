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
    
    % act based on data spacing
    if(leven(i)==h(v).true)
        data(i).x=detrend(data(i).x)+depmen(i);
    else
        data(i).amp=data(i).amp...
            -polyval(polyfit(data(i).times,data(i).amp,1),data(i).times)...
            +depmen(i);
    end
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end
