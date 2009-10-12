function [data]=removepolynomial(data,order)
%REMOVEPOLYNOMIAL    Remove polynomial trend from SEIZMO records
%
%    Usage:    data=removepolynomial(data,order)
%
%    Description: REMOVEPOLYNOMIAL(DATA,ORDER) removes the polynomial trend
%     of order ORDER from SEIZMO records.  For multi-component
%     records, each component is dealt with separately.  It is highly
%     recommended to combine this command with filtering operations to
%     reduce edge effects that may lead to poor data quality.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Check out the difference various order polynomials make:
%      plot1(removepolynomial(data(ones(1,10)),1:10))
%
%    See also: REMOVEMEAN, REMOVETREND, GETPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS

%     Version History:
%        June 24, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:30 GMT

% todo:

% check input
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% number of records
nrecs=numel(data);

% check order
if(~isnumeric(order) || any(order~=fix(order)) ...
        || ~any(numel(order)==[1 nrecs]))
    error('seizmo:removepolynomial:badOrder',...
        'ORDER must be a scalar or an array of integers.');
end
if(isscalar(order))
    order(1:nrecs,1)=order;
end

% header info
[b,e,npts]=getheader(data,'b','e','npts');
leven=getlgc(data,'leven');

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
        time=linspace(b(i),e(i),npts(i)).';
        for j=1:size(data(i).dep,2)
            data(i).dep(:,j)=data(i).dep(:,j) ...
                -polyval(polyfit(time,data(i).dep(:,j),order(i)),time);
        end
    % unevenly spaced
    else
        for j=1:size(data(i).dep,2)
            data(i).dep(:,j)=data(i).dep(:,j)-polyval(...
                polyfit(double(data(i).ind),data(i).dep(:,j),order(i)),...
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
