function [data]=ifourier(data)
%IFOURIER    Converts spectral SAClab records to the time domain
%
%    Description: Converts SAClab records from the frequency domain to the
%     time domain using the inverse fast fourier transform.  Output file
%     type is 'itime' (timeseries file).
%
%    Usage: data=ifourier(data)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: fourier, selctcmp

% check nargin
error(nargchk(1,1,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% retreive header info
[b,delta,sb,sdelta,nspts,leven,iftype]=gh(data,'b','delta','sb','sdelta','nspts','leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% loop through records
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of components
    ncmp=size(data(i).x);
    
    % check spacing and filetype
    if(leven(i)~=h(v).true)
        warning('SAClab:illegalOperation','Illegal operation on unevely spaced record %d',i);
        continue;
    elseif(iftype(i)==h(v).enum(1).val.ixyz)
        warning('SAClab:illegalOperation','Illegal operation on xyz file: record %d',i);
        continue;
    elseif(iftype(i)==h(v).enum(1).val.itime || iftype(i)==h(v).enum(1).val.ixy)
        warning('SAClab:illegalOperation','Illegal operation on timeseries/xy file: record %d',i);
        continue;
    elseif(isfield(h(v).enum(1).val,'incmp') && iftype(i)==h(v).enum(1).val.incmp)
       warning('SAClab:illegalOperation','Illegal operation on multi-component file: record %d',i);
        continue;
    elseif(iftype(i)~=h(v).enum(1).val.irlim && iftype(i)~=h(v).enum(1).val.iamph)
       warning( 'SAClab:illegalOperation','Illegal operation on record %d',i);
       continue;
    end
    
    % turn back into time domain
    if(iftype(i)==h(v).enum(1).val.irlim)
        data(i).x=nspts(i)*real(ifft(complex(data(i).x(:,1),data(i).x(:,2))));
        data(i)=ch(data(i),'iftype','itime');
    else
        data(i).x=nspts(i)/2*real(ifft(data(i).x(:,1).*exp(j*data(i).x(:,2))));
        data(i)=ch(data(i),'iftype','itime');
    end
    
    % truncate to original length and change class back
    data(i).x=oclass(data(i).x(1:nspts(i),:));
    
    % update header (note there is no field 'se')
    data(i)=ch(data(i),'b',sb(i),'e',sb(i)+(nspts(i)-1)*sdelta(i),...
        'delta',sdelta(i),'sb',b(i),'sdelta',delta(i),'nspts',ncmp(1),...
        'npts',nspts(i));
end

end