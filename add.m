function [data]=add(data,constant,cmp)
%ADD    Adds a constant to SAClab data records
%
%    Usage: [data]=add(data,constant,cmp)
%
%          - a scalar constant adds the value to all records
%          - a vector of constants (length must be same as number of 
%            records) allows adding different values to each record
%          - 'cmp' selects the dependent component to work on
%            (default is 1)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: subtr, mult, dvide

% check nargin
error(nargchk(1,3,nargin))

% no constant case
if(nargin==1 || isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% number of records
nrecs=length(data);

% check constant size
if(~isvector(constant) || ~any(length(constant)==[1 nrecs]))
    error('constant must be a scalar or vector of length equal to number of records')
end

% expand scalar constant
if(length(constant)==1); constant=constant(ones(nrecs,1)); end

% add constant and update header
for i=1:nrecs
    oclass=str2func(class(data(i).x));
    data(i).x(:,cmp)=oclass(double(data(i).x(:,cmp))+constant(i));
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
