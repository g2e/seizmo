function [data]=seizmofun(data,fun)
%SEIZMOFUN    Apply function to SEIZMO records
%
%    Description: SEIZMOFUN(DATA,FUN) applies the function defined by the 
%     function handle FUN to the dependent component(s) of SEIZMO records
%     in DATA.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record.
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, NCMP
%
%    Usage:    data=seizmofun(data,fun)
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 07:30 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seizmocheck(data,'dep'))

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

end
