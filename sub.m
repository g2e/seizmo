function [data]=sub(data,constant,cmp)
%SUB    Subtract a constant from SAClab data records
%
%    Description: SUB(DATA,CONSTANT) subtracts a constant from the 
%     dependent component(s) of SAClab data records.  For multi-component 
%     files, this operation is performed on every dependent component (this
%     includes spectral files).
%
%     SUB(DATA,CONSTANT,CMP) allows for operations on just components in
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
%    Usage:    data=sub(data,constant)
%              data=sub(data,constant,cmp_list)
%
%    Examples:
%     Do a Hilbert transform by converting records to the frequency 
%     domain, subtracting pi/2 from the phase (component 2 in spectral
%     records), and converting back to the time domain:
%      data=idft(sub(dft(data),pi/2,2))
%
%    See also: add, mul, divide, seisfun

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - improved input checks and docs
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        May  12, 2998 - dep* fix
%        June 12, 2008 - doc update, now works on all components by default
%        July 17, 2008 - history update, errors fixed, updated empty 
%                        component list behavior, .dep rather than .x, 
%                        allow constant to be an array, dataless support
%                        added, cmp checks, no longer uses add, doc update
%        Oct.  6, 2008 - minor code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  6, 2008 at 22:45 GMT

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
    error('SAClab:sub:badInput','Component list is bad!');
end

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('SAClab:sub:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('SAClab:sub:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% subtract constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    if(isempty(data(i).dep)); continue; end
    oclass=str2func(class(data(i).dep));
    data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))-constant(i));
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
