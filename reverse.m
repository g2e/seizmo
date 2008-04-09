function [data]=reverse(data)
%REVERSE    Time reverse SAClab data records
%
%    Description: Time reverses SAClab data records so that the beginning
%     and end of each record is switched.
%
%    Usage: [data]=reverse(data)
%
%    Examples:
%
%    See also: 

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% reverse records
for i=1:length(data)
    data(i).x=data(i).x(end:-1:1,:);
end

end
