function [data]=rmean(data)
%RMEAN    Remove data mean from SAClab data records
%
%    Description: RMEAN(DATA) removes the mean from data records.  This is
%     particularly useful for getting cleaner results from spectral
%     operations.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  data=rmean(data)
%
%    Examples:
%     It is generally a good idea to remove the mean from records before
%     performing any filtering operations like DECI to avoid edge effects:
%      p1(deci(data,5))        % more ringing
%      p1(deci(rmean(data),5)) % less ringing
%
%    See also: rtrend, taper

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Cleaned up documentation and added example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 02:55 GMT

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% remove mean and update header
for i=1:length(data)
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % loop through components
    for j=1:size(data(i).x,2)
        data(i).x(:,j)=data(i).x(:,j)-mean(data(i).x(:,j));
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % adjust header
    data(i)=ch(data(i),'depmen',mean(data(i).x(:)),...
        'depmin',min(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
