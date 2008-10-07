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
%        Oct.  3, 2008 - .dep & .ind
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct  3, 2008 at 15:15 GMT

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% remove mean and update header
for i=1:numel(data)
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % loop through components
    for j=1:size(data(i).dep,2)
        data(i).dep(:,j)=data(i).dep(:,j)-mean(data(i).dep(:,j));
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % adjust header
    data(i)=ch(data(i),'depmen',mean(data(i).dep(:)),...
        'depmin',min(data(i).dep(:)),'depmax',max(data(i).dep(:)));
end

end
