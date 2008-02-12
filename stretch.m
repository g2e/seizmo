function [data]=stretch(data,factor)
%STRETCH    Upsample SAClab data records
%
%    Description: Upsamples SAClab data records by an integer factor.  
%     Anti-aliasing is not an issue so there is no limit imposed on the
%     stretch factor.  Cascades are allowed regardless.  Uses the Matlab
%     function interp (Signal Processing Toolbox).
%
%    Usage:  [data]=stretch(data,factor)
%
%    Examples: 
%       to double the samplerates
%       [data]=stretch(data,2)
%
%       to cascade to a samplerate 40 times higher
%       [data]=stretch(data,[5 8])
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: deci, syncsr, intrplt
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
[delta,iftype,leven]=gh(data,'delta','iftype','leven');

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
        save=zeros((ncmp(1)-1)*of+1,ncmp(2));
        for j=1:ncmp(2)
            temp=data(i).x(:,j);
            % loop over factors
            for k=1:nf
                % stretch (filter length 4, cutoff freq is nyquist)
                temp=interp(temp,factor(k),4,1);
            end
            save(:,j)=temp;
        end
        
        % change class back
        data(i).x=oclass(save);
        
        % update header
        data(i)=ch(data(i),'npts',size(data(i).x,1),'delta',delta(i)/of,...
            'depmen',norm(mean(data(i).x)),'depmin',-norm(min(data(i).x)),...
            'depmax',norm(max(data(i).x)));
    end
end

end