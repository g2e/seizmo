function [data]=divide(data,constant,cmp)
%DIVIDE    Divide SEIZMO records by a constant
%
%    Usage:    data=divide(data,constant)
%              data=divide(data,constant,cmp_list)
%
%    Description:
%     DATA=DIVIDE(DATA,CONSTANT) divides the dependent component(s) of
%     SEIZMO records in DATA by a constant.  For multi-component files,
%     this operation is performed on each dependent component (including
%     spectral).
%
%     DATA=DIVIDE(DATA,CONSTANT,CMP) allows for operation on components in
%     the list CMP.  By default all components are operated on (use ':' to
%     replicate the default behavior).
%
%    Notes:
%     - A scalar constant applies the value to all records.
%     - A vector of constants (length must equal the number of records)
%       allows applying different values to each record.
%     - CMP is the dependent component(s) to work on (default is all).
%     - An empty list of components will not modify any components.
%    
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Alter the amplitudes of amplitude-phase spectral records without
%     % affecting the phase component by dividing only the first component:
%     data=divide(data,32,1)
%
%    See also: MULTIPLY, ADD, SUBTRACT, SOLOFUN

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
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Jan. 26, 2010 - seizmoverbose support, properly handle states
%        Jan.  6, 2011 - nargchk fix, seizmofun/solofun rename
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
    error('seizmo:divide:badInput','Component list is bad!');
end

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('seizmo:divide:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('seizmo:divide:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% detail message
if(verbose)
    disp('Dividing Record(s) by Constant')
    print_time_left(0,nrecs);
end

% divide by constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % divide selected cmp
    if(~isempty(cmp))
        oclass=str2func(class(data(i).dep));
        data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))/constant(i));
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
