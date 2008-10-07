function [data]=divide(data,constant,cmp)
%DIVIDE    Divide SAClab data records by a constant
%
%    Description: DIVIDE(DATA,CONSTANT) divides the dependent component(s)
%     of SAClab data records by a constant.  For multi-component files, 
%     this operation is performed on every dependent component (this 
%     includes spectral files).
%
%     DIVIDE(DATA,CONSTANT,CMP) allows for operation on just components in
%     the list CMP.  By default all components are operated on (use ':' to
%     replicate the default behavior).  See the examples section for a 
%     usage case.
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - CMP is the dependent component(s) to work on (default is all)
%     - an empty list of components will not modify any components
%
%    System requirements: Matlab 7
%    
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=divide(data,constant)
%              data=divide(data,constant,cmp_list)
%
%    Examples:
%     Alter the amplitudes of amplitude-phase spectral records without
%     affecting the phase component by dividing only the first component:
%      data=divide(data,32,1)
%
%    See also: mul, add, sub, seisfun

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - improved input checks and docs
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - doc update, now works on all components
%                        by default
%        June 20, 2008 - more doc updates
%        July  7, 2008 - history update, errors fixed, updated empty 
%                        component list behavior, .dep rather than .x, 
%                        allow constant to be an array, no longer uses mul
%        July 17, 2008 - doc update, dataless support added and cmp checks
%        Oct.  6, 2008 - minor code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2008 at 22:45 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin))

% check data structure
error(seischk(data,'dep'))

% no constant case
if(isempty(constant) || (nargin==3 && isempty(cmp))); return; end

% default component
if(nargin==2); cmp=':'; 
elseif(any(fix(cmp)~=cmp) || (~isnumeric(cmp) && ~strcmpi(':',cmp)))
    error('SAClab:divide:badInput','Component list is bad!');
end

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('SAClab:divide:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('SAClab:divide:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% divide by constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    oclass=str2func(class(data(i).dep));
    data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))/constant(i));
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
