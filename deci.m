function [data]=deci(data,factor)
%DECI    Downsample SAClab data records
%
%    Description: Downsamples SAClab data records by an integer factor.  
%     Keep the factor below 13 to avoid anti-aliasing filter problems.
%     Cascades are allowed by supplying multiple factors.  Uses the Matlab
%     function decimate (Signal Processing Toolbox).
%
%    Usage:  [data]=stretch(data,factors)
%
%    Examples: 
%     to halve the samplerates
%     [data]=deci(data,2)
%
%     to cascade to a samplerate 40 times lower
%     [data]=deci(data,[5 8])
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: stretch, syncsr, intrplt
%

% check inputs
error(nargchk(2,2,nargin))

% empty factor
if(isempty(factor)); return; end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% check factor
if(~isvector(factor) || factor<1)
    error('factor must be a vector of positive integer(s)')
end

% round factors
factor=round(factor);

% number of factors
nf=length(factor);

% overall factor
of=prod(factor);

% header info
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end
[b,delta,iftype,leven]=gh(data,'b','delta','iftype','leven');

% decimate and update header
for i=1:length(data)
    % header version
    v=data(i).version==vers;
    
    % must be evenly spaced
    if(leven(i)~=h(v).true)
        error('SAClab:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    elseif(iftype(i)~=h(v).enum(1).val.itime)
        error('SAClab:timeseriesOnly',...
            'Illegal operation on non-timeseries data');
    else
        % save class and convert to double precision
        oclass=str2func(class(data(i).x));
        data(i).x=double(data(i).x);
        
        % loop over components
        ncmp=size(data(i).x);
        for j=1:ncmp(2)
            temp=data(i).x(:,j);
            % loop over factors
            for k=1:nf
                % stretch (filter length 4, cutoff freq is nyquist)
                temp=decimate(temp,factor(k));
            end
            % preallocate
            if(j==1)
                save=zeros(size(temp,1),ncmp);
            end
            save(:,j)=temp;
        end
        
        % change class back
        data(i).x=oclass(save);
        
        % update header
        data(i)=ch(data(i),'npts',size(data(i).x,1),'delta',delta(i)*of,...
            'e',b(i)+(size(data(i).x,1)-1)*(delta(i)*of),...
            'depmen',norm(mean(data(i).x)),'depmin',-norm(min(data(i).x)),...
            'depmax',norm(max(data(i).x)));
    end
end

end