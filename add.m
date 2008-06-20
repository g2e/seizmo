function [data]=add(data,constant,cmp)
%ADD    Add a constant to SAClab data records
%
%    Description: ADD(DATA,CONSTANT) adds a constant to SAClab data
%     records.  For multi-component files, this operation is performed on
%     every component (this includes spectral files).
%
%     ADD(DATA,CONSTANT,CMP) allows for operations on just components in
%     the list CMP.  See the examples section for a usage case.
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - cmp_list gives the dependent component(s) to work on (default=all)
%
%    System requirements: Matlab 7
%
%    Data requirements: NONE
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=add(data,constant)
%           data=add(data,constant,cmp_list)
%
%    Examples:
%     Add a 135 degree (3*pi/4) phase shift to data records by only adding
%     to the phase component in amplitude-phase records (component 2):
%      data=idft(add(dft(data),3*pi/4,2))
%
%    See also: sub, mul, divide

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Cleaned up documentation and added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 00:20 GMT

% check nargin
error(nargchk(2,3,nargin))

% check data structure
error(seischk(data,'x'))

% no constant case
if(isempty(constant)); return; end

% default component
if(nargin==2 || isempty(cmp)); cmp=':'; end

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
    data(i)=ch(data(i),'depmax',max(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmen',mean(data(i).x(:)));
end

end
