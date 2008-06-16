function [data]=rtrend(data)
%RTREND    Remove linear trend from SAClab data records
%
%    Description: RTREND(DATA) removes the linear trend from SAClab data 
%     records by subtracting the best straight line fit to the data as 
%     determined by a least squares inversion.  It is highly recommended to
%     combine this command with any filtering operations to reduce edge
%     effects that may lead to poor data quality.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%    
%    Usage:  [data]=rtrend(data)
%
%    Examples:
%     4th order lowpass butter filter with a passband corner of 10s
%      data=iirfilter(rtrend(data),'low','butter',1/10,4)
%
%    See also: rmean, taper

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Cleaned up documentation and added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 03:15 GMT

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
