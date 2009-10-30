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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2009 at 17:00 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% try applying function
try
    % check input fun is a function
    if(~isa(fun,'function_handle'))
        error('seizmo:seizmofun:badInput','FUN must be a function handle!');
    end
    
    % apply function to records
    ncmp=nan(nrecs,1); npts=ncmp; depmen=ncmp; depmin=ncmp; depmax=ncmp;
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
    end
    
    % update header
    data=changeheader(data,'npts',npts,'ncmp',ncmp,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
