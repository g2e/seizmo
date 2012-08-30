function [data]=add(data,constant,cmp)
%ADD    Add a constant to SEIZMO records
%
%    Usage:    data=add(data,constant)
%              data=add(data,constant,cmp_list)
%
%    Description:
%     DATA=ADD(DATA,CONSTANT) adds a constant to the dependent component(s)
%     of SEIZMO records in DATA.  For multi-component files, this operation
%     is performed on every dependent component (including spectral).
%
%     DATA=ADD(DATA,CONSTANT,CMP) allows for operation on components in the
%     list CMP.  By default all components are operated on (use ':' to
%     replicate the default behavior).
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
%     % Add a 135 degree (3*pi/4) phase shift to records by only adding
%     % to the phase component in amplitude-phase records (component 2):
%     data=idft(add(dft(data),3*pi/4,2))
%
%    See also: SUBTRACT, MULTIPLY, DIVIDE, SOLOFUN

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
%        Jan. 26, 2010 - seizmoverbose support, properly handle states
%        Aug. 16, 2010 - nargchk fix, versioninfo instead of seizmocheck
%        Jan.  6, 2011 - back to seizmocheck, seizmofun/solofun rename
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 18:25 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% no constant case
if(isempty(constant) || (nargin==3 && isempty(cmp))); return; end

% default component
if(nargin==2); cmp=':'; 
elseif(any(fix(cmp)~=cmp) || (~isnumeric(cmp) && ~strcmpi(':',cmp)))
    error('seizmo:add:badInput','Component list is bad!');
end

% verbosity
verbose=seizmoverbose;

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

% detail message
if(verbose)
    disp('Adding Constant to Record(s)');
    print_time_left(0,nrecs);
end

% add constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % add to selected cmp
    if(~isempty(cmp))
        oclass=str2func(class(data(i).dep));
        data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))+constant(i));
    end
    
    % dep*
    depmen(i)=nanmean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% try updating header
try
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
