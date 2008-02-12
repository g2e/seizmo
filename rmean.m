function [data]=rmean(data)
%RMEAN    Remove data mean from SAClab data records
%
%    Usage:  [data]=rmean(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: rslope, rdrift

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

% remove mean and update header
for i=1:length(data)
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % loop through components
    for j=1:size(data(i).x,2)
        data(i).x(:,j)=data(i).x(:,j)-mean(data(i).x(:,j));
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end
