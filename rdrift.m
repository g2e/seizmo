function [data]=rdrift(data)
%RDRIFT    Remove mean and linear trend from seislab data records
%
%    Description: Removes the mean and trend from records.  Uses the Matlab
%     functions detrend or polyfit and polyval.
%
%    Usage:  [data]=rdrift(data)
%
%    See also: rmean, rslope

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% header info
leven=glgc(data,'leven');

% check leven
t=strcmp(leven,'true');
f=strcmp(leven,'false');
if(~all(t | f))
    error('sieslab:rdrift:levenBad',...
        'logical field leven needs to be set'); 
end

% remove trend and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % act based on data spacing
    if(strcmp(leven(i),'true'))
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=detrend(data(i).x(:,j));
        end
    else
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=data(i).x(:,j) ...
                -polyval(polyfit(double(data(i).t),data(i).x(:,j),1),...
                double(data(i).t));
        end
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)));
end

end
