function [data]=seizmofun(data,fun)
%SEIZMOFUN    Apply function to SEIZMO records
%
%    Usage:    data=seizmofun(data,fun)
%
%    Description: SEIZMOFUN(DATA,FUN) applies the function defined by the 
%     function handle FUN to the dependent component(s) of SEIZMO records
%     in DATA.  Does not adjust the independent data at all.
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
%     A multi-tool of a function:
%      data=seizmofun(data,@abs)
%      data=seizmofun(data,@sign)
%      data=seizmofun(data,@log)
%      data=seizmofun(data,@sqrt)
%      data=seizmofun(data,@exp)
%      data=seizmofun(data,@(x)log(x)/log(4))
%      data=seizmofun(data,@(x)x.^3)
%      data=seizmofun(data,@(x)3.^x)
%      data=seizmofun(data,@(x)real(exp(-2*i*pi*x)))
%
%    See also: ADD, DIVIDE, MULTIPLY, SUBTRACT, SLIDINGFUN

%     Version History:
%        Apr.  9, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 17, 2008 - doc update, dataless support, .dep rather than .x,
%                        added history, single ch call
%        Oct.  6, 2008 - use CHKHDR over CH (slims code)
%        Nov. 22, 2008 - update for new name schema (now SEIZMOFUN)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct. 15, 2009 - no header checking, updates header for npts, ncmp,
%                        dep*, catches errors to keep checking correct
%        Dec.  4, 2009 - whoops, forgot NRECS
%        Feb.  3, 2010 - proper SEIZMO handling, versioninfo caching,
%                        seizmoverbose support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2010 at 17:05 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% try applying function
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % check input fun is a function
    if(~isa(fun,'function_handle'))
        error('seizmo:seizmofun:badInput',...
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
        oclass=str2func(class(data(i).dep));
        data(i).dep=oclass(fun(double(data(i).dep)));
        
        % get npts, ncmp, dep*
        [npts(i),ncmp(i)]=size(data(i).dep);
        if(npts(i)) % skip dep* for empty
            depmen(i)=mean(data(i).dep(:)); 
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
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
