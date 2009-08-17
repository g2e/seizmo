function [data]=multiply(data,constant,cmp)
%MULTIPLY    Multiply SEIZMO records by a constant
%
%    Usage:    data=multiply(data,constant)
%              data=multiply(data,constant,cmp_list)
%
%    Description: MULTIPLY(DATA,CONSTANT) multiplies the dependent
%     component(s) of SEIZMO records by a constant.  For multi-component
%     files, this operation is performed on every dependent component (this
%     includes spectral files).
%
%     MULTIPLY(DATA,CONSTANT,CMP) allows for operations on just components
%     in the list CMP.  By default all components are operated on (use ':'
%     to replicate the default behavior).  See the examples section for a 
%     usage case.
%
%    Notes:
%     - a scalar constant applies the value to all records
%     - a vector of constants (length must equal the number of records)
%       allows applying different values to each record
%     - CMP is the dependent component(s) to work on (default is all)
%     - an empty list of components will not modify any components
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Get the complex conjugate of a real-imaginary spectral records by
%     multiplying the imaginary component by -1 (component 2):
%      data=multiply(data,-1,2)
%
%    See also: subtract, add, divide, seizmofun

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - improved input checks and docs
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - doc update, now works on all components by default
%        July 17, 2008 - history update, errors fixed, updated empty 
%                        component list behavior, .dep rather than .x, 
%                        allow constant to be an array, dataless support,
%                        cmp checks, and doc update
%        Oct.  6, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema (now MULTIPLY)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:30 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% no constant case
if(isempty(constant) || (nargin==3 && isempty(cmp))); return; end

% default component
if(nargin==2); cmp=':'; 
elseif(any(fix(cmp)~=cmp) || (~isnumeric(cmp) && ~strcmpi(':',cmp)))
    error('seizmo:multiply:badInput','Component list is bad!');
end

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('seizmo:multiply:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('seizmo:multiply:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% multiply by constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    if(~isempty(cmp))
        oclass=str2func(class(data(i).dep));
        data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))*constant(i));
    end
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
