function [data]=dif(data)
%DIF    Differentiates SAClab data records using discrete differences
%
%    Description: Calculates and returns the derivative of each record
%     using the differences between points as an approximation of the 
%     derivative at the midpoint.  Works with unevenly spaced data.  Uses
%     the Matlab function diff.
%
%    Notes: 
%       Timing is shifted to midpoints.
%       Reduces npts by one.
%
%    Usage: [data]=dif(data);
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: integrt, integrt2

% check nargin
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

% retreive header info
[delta,b,e,npts,leven]=gh(data,'delta','b','e','npts','leven');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% take derivative and update header
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % evenly spaced?
    if (leven(i)==h(v).true)
        data(i).x=diff(data(i).x)/delta(i);
        b(i)=b(i)+delta(i)/2; e(i)=e(i)-delta(i)/2; npts(i)=npts(i)-1;
    else
        data(i).t=double(data(i).t);
        t=diff(data(i).t);
        data(i).x=diff(data(i).x)./t(:,ones(1,size(data(i).x,2)));
        data(i).t=oclass(data(i).t(1:npts-1)+t/2);
        npts(i)=npts(i)-1; b(i)=data(i).t(1); e(i)=data(i).t(end);
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)),...
        'b',b(i),'e',e(i),'npts',npts(i));
end

end