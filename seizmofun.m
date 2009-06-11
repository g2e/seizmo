function [data]=seizmofun(data,fun)
%SEIZMOFUN    Apply function to SEIZMO records
%
%    Usage:    data=seizmofun(data,fun)
%
%    Description: SEIZMOFUN(DATA,FUN) applies the function defined by the 
%     function handle FUN to the dependent component(s) of SEIZMO records
%     in DATA.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record.
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
%    See also: add, divide, multiply, subtract, slidingfun

%     Version History:
%        Apr.  9, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 17, 2008 - doc update, dataless support, .dep rather than .x,
%                        added history, single ch call
%        Oct.  6, 2008 - use CHKHDR over CH (slims code)
%        Nov. 22, 2008 - update for new name schema (now SEIZMOFUN)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June  9, 2009 - added testing table
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
%     Last Updated June  9, 2009 at 19:35 GMT

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

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('seizmo:seizmofun:badInput','FUN must be a function handle!');
end

% apply function to records
for i=1:numel(data)
    oclass=str2func(class(data(i).dep));
    data(i).dep=oclass(fun(double(data(i).dep)));
end

% update header
data=checkheader(data);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
