function [data]=solofun(data,fun)
%SOLOFUN    Apply function to SEIZMO records individually
%
%    Usage:    data=solofun(data,fun)
%
%    Description:
%     SOLOFUN(DATA,FUN) applies the function defined by the function handle
%     FUN to the dependent component(s) of SEIZMO records in DATA.  Does
%     not adjust the independent data at all.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record.
%     - Changing the NPTS in an unevenly sampled record will not adjust the
%       independent component data.  You must do this yourself!
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, NCMP
%
%    Examples:
%     % A multi-tool of a function:
%     data=solofun(data,@abs)
%     data=solofun(data,@sign)
%     data=solofun(data,@log)
%     data=solofun(data,@sqrt)
%     data=solofun(data,@exp)
%     data=solofun(data,@(x)log(x)/log(4))
%     data=solofun(data,@(x)x.^3)
%     data=solofun(data,@(x)3.^x)
%     data=solofun(data,@(x)real(exp(-2*i*pi*x)))
%
%    See also: ADD, DIVIDE, MULTIPLY, SUBTRACT, SLIDINGFUN

%     Version History:
%        Apr.  9, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 17, 2008 - doc update, dataless support, .dep rather than .x,
%                        added history, single ch call
%        Oct.  6, 2008 - use CHKHDR over CH (slims code)
%        Nov. 22, 2008 - update for new name schema (now SOLOFUN)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct. 15, 2009 - no header checking, updates header for npts, ncmp,
%                        dep*, catches errors to keep checking correct
%        Dec.  4, 2009 - whoops, forgot NRECS
%        Feb.  3, 2010 - proper SEIZMO handling, versioninfo caching,
%                        seizmoverbose support
%        Jan.  6, 2011 - renamed from seizmofun to solofun to better
%                        distinguish from multifun (was recordfun), removed
%                        versioninfo caching (too prone to breakage)
%        Feb.  5, 2012 - doc update, code comment
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2012 at 17:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% try applying function
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check input fun is a function
    if(~isa(fun,'function_handle'))
        error('seizmo:solofun:badInput',...
            'FUN must be a function handle!');
    end
    
    % detail message
    if(verbose)
        disp('Applying Function to the Dependent Data of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % apply function to records
    ncmp=nan(nrecs,1); npts=ncmp;
    depmen=ncmp; depmin=ncmp; depmax=ncmp;
    for i=1:nrecs
        % save class, convert to double precision,
        % apply function & convert back to original precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=oclass(fun(double(data(i).dep)));
        
        % get npts, ncmp, dep*
        [npts(i),ncmp(i)]=size(data(i).dep);
        if(npts(i)) % skip dep* for empty
            depmen(i)=nanmean(data(i).dep(:)); 
            depmin(i)=min(data(i).dep(:)); 
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update header
    data=changeheader(data,'npts',npts,'ncmp',ncmp,...
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
