function [data]=add(data,constant,cmp)
%ADD    Add a constant to SEIZMO records
%
%    Usage:    data=add(data,constant)
%              data=add(data,constant,cmp_list)
%
%    Description: ADD(DATA,CONSTANT) adds a constant to the dependent 
%     component(s) of SEIZMO records.  For multi-component files, this
%     operation is performed on every dependent component (this includes 
%     spectral files).
%
%     ADD(DATA,CONSTANT,CMP) allows for operation on just components in
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
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Add a 135 degree (3*pi/4) phase shift to records by only adding
%     to the phase component in amplitude-phase records (component 2):
%      data=idft(add(dft(data),3*pi/4,2))
%
%    See also: subtract, multiply, divide, seizmofun

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - improved input checks and docs
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - doc update, now works on all components
%                        by default
%        June 20, 2008 - more doc updates
%        June 28, 2008 - history update, ch only called once now, errors
%                        fixed, updated empty component list behavior,
%                        .dep rather than .x
%        July  7, 2008 - allow constant to be an array
%        July 17, 2008 - doc update, dataless support added and cmp checks
%        Oct.  6, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schemas
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 29, 2009 - add testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 29, 2009 at 01:30 GMT

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
    error('seizmo:add:badInput','Component list is bad!');
end

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('seizmo:add:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('seizmo:add:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% add constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    if(~isempty(cmp))
        oclass=str2func(class(data(i).dep));
        data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))+constant(i));
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
