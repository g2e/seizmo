function [data]=removetrend(data)
%REMOVETREND    Remove linear trend from SEIZMO records
%
%    Description: REMOVETREND(DATA) removes the linear trend from SEIZMO
%     records by subtracting the best straight line fit to the data as 
%     determined by a least squares inversion.  For multi-component
%     records, each component is dealt with separately.  It is highly
%     recommended to combine this command with any filtering operations to
%     reduce edge effects that may lead to poor data quality.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%    
%    Usage:    data=removetrend(data)
%
%    Examples:
%     4th order lowpass butter filter with a passband corner of 10s
%      data=iirfilter(removetrend(data),'low','butter',1/10,4)
%
%    See also: removemean, taper, removedeadrecords

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Nov. 27, 2007 - RDRIFT added (removes both mean and trend)
%        Jan. 14, 2008 - handle uneven files
%        Feb. 12, 2008 - RTREND renamed to RSLOPE
%        Feb. 29, 2008 - SEISCHK support, handle uneven files better,
%                        RTREND removed (RSLOPE still around)
%        Mar.  4, 2008 - doc update, minor code cleaning
%        May  12, 2008 - fix dep* formula
%        June 12, 2008 - doc update, history added, renamed from RDRIFT to
%                        RTREND to match SAC, dropped RSLOPE
%        Oct.  3, 2008 - .dep & .ind
%        Nov. 22, 2008 - doc update, history fix, renamed from RTREND to
%                        REMOVETREND, one CHANGEHEADER call, better checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 09:05 GMT

% todo:

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% header info
leven=getlgc(data,'leven');

% number of records
nrecs=numel(data);

% remove trend and update header
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:numel(data)
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
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
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% adjust header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
