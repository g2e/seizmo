function [data]=subtract(data,constant,cmp)
%SUBTRACT    Subtract a constant from SEIZMO records
%
%    Usage:    data=subtract(data,constant)
%              data=subtract(data,constant,cmp_list)
%
%    Description:
%     DATA=SUBTRACT(DATA,CONSTANT) subtracts a constant from the dependent
%     component(s) of SEIZMO records in DATA.  For multi-component files,
%     this operation is performed on every dependent component (including
%     spectral).
%
%     DATA=SUBTRACT(DATA,CONSTANT,CMP) allows for operations on components
%     in the list CMP.  By default all components are operated on (use ':'
%     to replicate the default behavior).
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
%     % Do a Hilbert transform by converting records to the frequency 
%     % domain, subtracting pi/2 from the phase (component 2 in spectral
%     % records), and converting back to the time domain:
%     data=idft(subtract(dft(data),pi/2,2))
%
%    See also: ADD, MULTIPLY, DIVIDE, SOLOFUN

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
%        Nov. 22, 2008 - update to new name schema (now SUBTRACT)
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
    error('seizmo:subtract:badInput','Component list is bad!');
end

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% check constant
if(~isnumeric(constant))
    error('seizmo:subtract:badInput','Constant must be numeric!');
elseif(isscalar(constant))
    constant=constant(ones(nrecs,1));
elseif(numel(constant)~=nrecs)
    error('seizmo:subtract:badInput',...
        'Number of elements in constant not equal to number of records!');
end

% detail message
if(verbose)
    disp('Subtracting Constant from Record(s)')
    print_time_left(0,nrecs);
end

% subtract constant
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % subtract from selected cmp
    if(~isempty(cmp))
        oclass=str2func(class(data(i).dep));
        data(i).dep(:,cmp)=oclass(double(data(i).dep(:,cmp))-constant(i));
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
