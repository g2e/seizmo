function [data]=integrt(data)
%INTEGRT    Integrates SAClab data records using the trapezoidal rule
%
%    Description: Calculates and returns the integral of each record using
%     the trapezoidal rule.  Works with unevenly spaced data.  Uses the 
%     Matlab function cumtrapz.
%
%    Notes:
%       No change to timing. 
%       No change to npts.
%       First point is 0.
%
%    Usage: [data]=integrt(data);
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: dif, integrt2

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
[delta,leven]=gh(data,'delta','leven');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% integrate and update header
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % evenly spaced?
    if (leven(i)==h(v).true)
        data(i).x=delta(i)*cumtrapz(data(i).x);
    else
        data(i).x=cumtrapz(double(data(i).t),data(i).x);
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end