function [p]=getpolynomial(data,order)
%GETPOLYNOMIAL    Get polynomial fit to SEIZMO records
%
%    Usage:    p=getpolynomial(data,order)
%
%    Description: P=GETPOLYNOMIAL(DATA,ORDER) gets the polynomial trend
%     of order ORDER from SEIZMO records.  For multi-component records,
%     each component is dealt with separately.  Polynomial coefficients are
%     returned in cell array P, with each cell corresponding to each record
%     in DATA and each row in each cell giving the coefficients for each
%     component.
%
%    Notes:
%
%    Examples:
%     Get various polynomial fits to a record:
%      getpolynomial(data(ones(1,5)),0:4)
%
%    See also: REMOVEMEAN, REMOVETREND, REMOVEPOLYNOMIAL, TAPER,
%              REMOVEDEADRECORDS

%     Version History:
%        June 24, 2009 - initial version
%        Dec.  4, 2009 - drop linspace usage (speed/accuracy decision)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  4, 2009 at 18:55 GMT

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
[delta,npts]=getheader(data,'delta','npts');
leven=getlgc(data,'leven');

% remove trend and update header
p=cell(nrecs,1);
for i=1:numel(data)
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % convert to double precision
    data(i).dep=double(data(i).dep);
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
        time=((0:npts(i)-1)*delta(i)).';
        for j=1:size(data(i).dep,2)
            p{i}(j,1:order(i)+1)=polyfit(time,data(i).dep(:,j),order(i));
        end
    % unevenly spaced
    else
        for j=1:size(data(i).dep,2)
            p{i}(j,1:order(i)+1)=...
                polyfit(double(data(i).ind),data(i).dep(:,j),order(i));
        end
    end
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
