function [data]=integrt(data)
%INTEGRT    Integrates SAClab data records using the trapezoidal rule
%
%    Description: Calculates and returns the integral of SAClab data
%     records using the trapezoidal rule.  Works with unevenly spaced data.
%     Uses the Matlab function cumtrapz.
%
%    Notes:
%     - No change to timing. 
%     - No change to npts.
%     - First point is 0.
%
%    Usage: [data]=integrt(data)
%
%    Examples:
%
%    See also: dif, integrt2

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% retreive header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
delta=gh(data,'delta');

% integrate and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % evenly spaced?
    if (strcmp(leven(i),'true'))
        data(i).x=delta(i)*cumtrapz(data(i).x);
    else
        data(i).x=cumtrapz(double(data(i).t),data(i).x);
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end
