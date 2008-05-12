function [data]=rdrift(data)
%RDRIFT    Remove mean and linear trend from SAClab data records
%
%    Description: Removes the mean and trend from SAClab data records.  
%     Uses Matlab functions detrend, polyfit and polyval.
%
%    Usage:  [data]=rdrift(data)
%
%    Examples:
%
%    See also: rmean, rslope

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))

% remove trend and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=detrend(data(i).x(:,j));
        end
    % unevenly spaced
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
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
