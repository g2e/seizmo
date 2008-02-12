function [data]=subtr(data,constant,cmp)
%SUBTR    Subtracts a constant to SAClab data records
%
%    Usage: [data]=subtr(data,constant,cmp)
%
%          - a scalar constant subtracts the value to all records
%          - a vector of constants (length must be same as number of 
%            records) allows subtracting different values to each record
%          - 'cmp' selects the dependent component to work on
%            (default is 1)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: add, mult, dvide

% check nargin
error(nargchk(1,3,nargin))

% no constant case
if(nargin==1 || isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% send to add
data=add(data,-constant,cmp);

end