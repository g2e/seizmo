function [data]=dif(data)
%DIF    Differentiate seislab data records using discrete differences
%
%    Description: Calculates and returns the derivative of each record
%     using the differences between points as an approximation of the 
%     derivative at the midpoint.  Works with unevenly spaced data.  Uses
%     the Matlab function diff.
%
%    Notes: 
%     - Timing is shifted to midpoints.
%     - Reduces npts by one.
%
%    Usage: [data]=dif(data);
%
%    See also: integrt, integrt2

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% retreive header info
leven=glgc(data,'leven');
[delta,b,e,npts]=gh(data,'delta','b','e','npts');

% take derivative and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % evenly spaced?
    if(strcmp(leven(i),'true'))
        data(i).x=diff(data(i).x)/delta(i);
        b(i)=b(i)+delta(i)/2; e(i)=e(i)-delta(i)/2; npts(i)=npts(i)-1;
    elseif(strcmp(leven(i),'false'))
        data(i).t=double(data(i).t);
        t=diff(data(i).t);
        data(i).x=diff(data(i).x)./t(:,ones(1,size(data(i).x,2)));
        data(i).t=oclass(data(i).t(1:npts-1)+t/2);
        npts(i)=npts(i)-1; b(i)=data(i).t(1); e(i)=data(i).t(end);
    else
        error('sieslab:dif:levenBad',...
            'logical field leven needs to be set for record %d',i);
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)),...
        'b',b(i),'e',e(i),'npts',npts(i));
end

end

