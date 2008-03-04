function [data]=add(data,constant,cmp)
%ADD    Add a constant to SAClab data records
%
%    Usage: [data]=add(data,constant,cmp)
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - 'cmp' selects the dependent component to work on (default is 1)
%
%    Examples:
%
%    See also: subtr, mult, dvide

% check nargin
error(nargchk(2,3,nargin))

% check data structure
error(seischk(data,'x'))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% number of records
nrecs=length(data);

% check constant
if(~isnumeric(constant)); error('constant must be numeric');
elseif(isscalar(constant)); constant=constant(ones(nrecs,1));
elseif(~isvector(constant) || length(constant)~=nrecs)
    error('constant vector length ~= # records')
end

% add constant and update header
for i=1:nrecs
    oclass=str2func(class(data(i).x));
    data(i).x(:,cmp)=oclass(double(data(i).x(:,cmp))+constant(i));
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
