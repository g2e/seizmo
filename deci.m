function [data]=deci(data,factor)
%DECI    Downsample SAClab data records by integer factor
%
%    Description: DECI(DATA,FACTOR) downsamples SAClab data records by an 
%     integer FACTOR after implementing an anti-aliasing filter.  To avoid
%     adding significant numerical noise to the data, keep the decimation
%     factor below 13.  Strong amplitudes at or near the start and end of 
%     records will introduce edge effects that can be reduced by first 
%     detrending and tapering records.  Cascades are allowed by supplying 
%     a vector of factors.  Uses the Matlab function decimate (Signal 
%     Processing Toolbox).
%
%    Header changes: DELTA, NPTS, E, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  data=deci(data,factors)
%
%    Examples: 
%     To halve the samplerates of records in data:
%      [data]=deci(data,2)
%
%     To cascade records to a samplerate 40 times lower;
%      [data]=deci(data,[5 8])
%
%    See also: stretch, syncsr, interpolate

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 03:45 GMT

% check inputs
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% empty factor
if(isempty(factor)); return; end

% check factor
if(~isnumeric(factor) || ~isvector(factor) || ...
        any(factor<1) || any(fix(factor)~=factor))
    error('SAClab:deci:badInput',...
        'factor must be a scalar/vector of positive integer(s)')
end

% number of factors
nf=length(factor);

% overall factor
of=prod(factor);

% check spacing
if(any(~strcmp(glgc(data,'leven'),'true')))
    error('SAClab:deci:evenlySpacedOnly',...
        'illegal operation on unevenly spaced data');
end

% header info
[b,delta]=gh(data,'b','delta');
delta=delta*of;

% decimate and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % loop over components
    ncmp=size(data(i).x,2);
    for j=1:ncmp
        temp=data(i).x(:,j);
        % loop over factors
        for k=1:nf
            temp=decimate(temp,factor(k));
        end
        % preallocate after we know length
        if(j==1)
            len=size(temp,1);
            save=zeros(len,ncmp);
        end
        save(:,j)=temp;
    end
    
    % change class back
    data(i).x=oclass(save);
    
    % update header
    data(i)=ch(data(i),'npts',len,'delta',delta(i),...
        'e',b(i)+(len-1)*delta(i),'depmin',min(data(i).x(:)),...
        'depmen',mean(data(i).x(:)),'depmax',max(data(i).x(:)));
end

end
