function [data]=mul(data,constant,cmp)
%MUL    Multiply SAClab data records by a constant
%
%    Description: MUL(DATA,CONSTANT) multiplies the dependent component(s) 
%     of SAClab data records by a constant.  For multi-component files, 
%     this operation is performed on every dependent component (this 
%     includes spectral files).
%
%     MUL(DATA,CONSTANT,CMP) allows for operations on just components in
%     the list CMP.  By default all components are operated on (use ':' to
%     replicate the default behavior).  See the examples section for a 
%     usage case.
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - cmp_list is the dependent component(s) to work on (default is all)
%     - an empty list of components will not modify any components
%
%    System requirements: Matlab 7
%
%    Data requirements: NONE
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=mul(data,constant)
%           data=mul(data,constant,cmp_list)
%
%    Examples:
%     Get the complex conjugate of a real-imaginary spectral records by
%     multiplying the imaginary component by -1 (component 2):
%      data=mul(data,-1,2)
%
%    See also: sub, add, divide, seisfun

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - improved input checks and documentation
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor documentation update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - documentation update, now works on all components
%                        by default
%        July 17, 2008 - history update, errors fixed, updated empty 
%                        component list behavior, .dep rather than .x, 
%                        allow constant to be an array, dataless support,
%                        cmp checks, and documentation update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2008 at 08:10 GMT

% todo:
%

% check nargin
error(nargchk(2,3,nargin))

% check data structure
error(seischk(data,'dep'))

% no constant case
if(isempty(constant) || (nargin==3 && isempty(cmp))); return; end

% default component
if(nargin==2); cmp=':'; 
elseif(any(fix(cmp)~=cmp) || (~isnumeric(cmp) && ~strcmpi(':',cmp)))
    error('SAClab:mul:badInput','Component list is bad');
end

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('SAClab:mul:badInput','constant must be numeric');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('SAClab:mul:badInput',...
        'number of elements in constant not equal to number of records')
end

% multiply by constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    oclass=str2func(class(data(i).dep));
    data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))*constant(i));
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
