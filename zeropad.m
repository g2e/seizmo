function [data]=zeropad(data,padpow)
%ZEROPAD    Extend Saclab data records to a 2^n length by zero padding
%
%    Description: Adds zeros to the end of records to extend their length
%     to a power of 2.  This is primarily for frequency domain work as it
%     makes the discrete fourier transform process faster.  Higher powers
%     give more frequency precision (but not necessarily better accuracy!).
%
%    Usage: [data]=zeropad(data,padpow)
%
%    Examples:
%     Pad records so that they are the next higher power of 2 in length
%      data=zeropad(data,0)
%     Pad records with zeros to double the length of the last example
%      data=zeropad(data,1)
%     Truncate records to next shorter power of 2 length
%      data=zeropad(data,-1)
%
%    See also: fourier, ifourier, iirfilter

% check nargin
error(nargchk(2,2,nargin))

% check data structure
if(~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% number of records
nrecs=length(data);

% retreive header info
[leven]=gh(data,'leven');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% no constant case
if(all(round(padpow)~=padpow) || isempty(padpow) || ...
        ~isvector(padpow) || ~any(length(padpow)==[1 nrecs]));
    error('padding argument must be an integer'); 
end

% expand scalar constant
if(length(padpow)==1); padpow=padpow(ones(nrecs,1)); end

% add zeros and update header
for i=1:nrecs
    % logical index of header
    v=data(i).version==vers;
    
    % skip unevenly spaced files
    if(leven(i)==h(v).false)
        warning('SAClab:zeropad:illegalOperation',...
            'Illegal operation with unevenly sampled record %g',i);
        continue;
    end
    
    % pad and update
    sx=size(data(i).x,1);
    newlen=2^(nextpow2(sx(1))+padpow(i));
    data(i).x=[data(i).x; zeros(newlen-sx(1),sx(2),class(data(i).x))];
    data(i)=ch(data(i),'depmax',norm(max(data(i).x)),'npts',newlen,...
        'depmin',-norm(min(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end