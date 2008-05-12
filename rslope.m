function [data]=rslope(data)
%RSLOPE    Remove slope from SAClab data records
%
%    Description: Removes the slope from records while preserving the data
%     mean.  Uses Matlab functions detrend, polyfit and polyval.
%
%    Usage:  [data]=rslope(data)
%
%    Examples:
%
%    See also: rmean, rdrift

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% grab timing info
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
            data(i).x(:,j)=detrend(data(i).x(:,j))+mean(data(i).x(:,j));
        end
    % unevenly spaced
    else
        for j=1:size(data(i).x,2)
            data(i).x(:,j)=data(i).x(:,j) ...
                -polyval(polyfit(double(data(i).t),data(i).x(:,j),1),...
                double(data(i).t))+mean(data(i).x(:,j));
        end
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
