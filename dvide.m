function [data]=dvide(data,constant,cmp)
%DVIDE    Divides SAClab data records by a constant
%
%    Usage: [data]=dvide(data,constant,cmp)
%
%          - a scalar constant applies the value to all records
%          - a vector of constants (length must be same as number of 
%            records) allows applying different values to each record
%          - 'cmp' selects the dependent component to work on
%            (default is 1)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: mult, add, subtr

% check nargin
error(nargchk(1,3,nargin))

% no constant case
if(nargin==1 || isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% send to mult
data=mult(data,1./constant,cmp);

end