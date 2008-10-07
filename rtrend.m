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
%        Oct.  3, 2008 - .dep & .ind
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct  3, 2008 at 15:15 GMT

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))

% remove trend and update header
for i=1:numel(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
        for j=1:size(data(i).dep,2)
            data(i).dep(:,j)=detrend(data(i).dep(:,j));
        end
    % unevenly spaced
    else
        for j=1:size(data(i).dep,2)
            data(i).dep(:,j)=data(i).dep(:,j) ...
                -polyval(polyfit(double(data(i).ind),data(i).dep(:,j),1),...
                double(data(i).ind));
        end
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % adjust header
    data(i)=ch(data(i),'depmen',mean(data(i).dep(:)),...
        'depmin',min(data(i).dep(:)),'depmax',max(data(i).dep(:)));
end

end
