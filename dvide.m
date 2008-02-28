function [data]=dvide(data,constant,cmp)
%DVIDE    Divide seislab data records by a constant
%
%    Usage: [data]=dvide(data,constant,cmp)
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - 'cmp' selects the dependent component to work on (default is 1)
%
%    See also: mult, add, subtr

% check nargin
error(nargchk(2,3,nargin))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=1; end

% send to mult
data=mult(data,1./constant,cmp);

end
