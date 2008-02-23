function [data]=subtr(data,constant,cmp)
%SUBTR    Subtract a constant from SAClab data records
%
%    Usage: [data]=subtr(data,constant,cmp)
%
%          - a scalar constant applies the value to all records
%          - a vector of constants (length must equal the number of 
%            records) allows applying different values to each record
%          - 'cmp' selects the dependent component to work on
%            (default is 1)
%
%    See also: add, mult, dvide

% check nargin
error(nargchk(2,3,nargin))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% send to add
data=add(data,-constant,cmp);

end
