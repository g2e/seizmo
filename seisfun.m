function [data]=seisfun(data,fun)
%SEISFUN    Apply function to SAClab data records
%
%    Description: SEISFUN(DATA,FUN) applies function FUN to SAClab data 
%     records in DATA.  Function FUN is expected to only modify the
%     amplitudes of the records and not the length or timing.  The number
%     of components in a multi-component record can be changed.
%
%    Usage: [data]=seisfun(data,fun)
%
%    Note: Will update header values depmen, depmin, depmax.  Does not
%     check the number of points in the resulting record and thus does not
%     update the npts or e header values.
%
%    Examples:
%
%     The multi-tool of functions:
%      seisfun(data,@abs)
%      seisfun(data,@sign)
%      seisfun(data,@log)
%      seisfun(data,@sqrt)
%      seisfun(data,@exp)
%      seisfun(data,@(x)log(x)/log(4))
%      seisfun(data,@(x)x.^3)
%      seisfun(data,@(x)3.^x)
%      seisfun(data,@(x)real(exp(-2*i*pi*x)))
%
%    See also: 

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('FUN must be a function handle!')
end

% apply function to records
for i=1:length(data)
    % preserve class while passing data to fun as double
    oclass=str2func(class(data(i).x));
    data(i).x=oclass(fun(double(data(i).x)));
    
    % update header
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),...
        'depmin',-norm(min(data(i).x)),...
        'depmen',norm(mean(data(i).x)));
end

end
