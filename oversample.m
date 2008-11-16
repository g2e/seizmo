function [data]=stretch(data,factor)
%STRETCH    Upsample SEIZMO data records by an integer factor
%
%    Description: STRETCH(DATA,FACTOR) upsamples SEIZMO data records by an 
%     integer FACTOR.  Anti-aliasing is not an issue so there is no limit 
%     imposed on the stretch factor.  Cascades are allowed regardless.  
%     Uses the Matlab function interp (Signal Processing Toolbox).
%
%    Header changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  [data]=stretch(data,factor)
%
%    Examples: 
%     To double samplerates:
%      data=stretch(data,2)
%
%     To cascade to a samplerate 40 times higher:
%      data=stretch(data,[5 8])
%
%    See also: deci, syncsr, interpolate

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 03:55 GMT

% check inputs
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% empty factor
if(isempty(factor)); return; end

% check factor
if(~isnumeric(factor) || ~isvector(factor) || ...
        any(factor<1) || any(fix(factor)~=factor))
    error('seizmo:stretch:badInput',...
        'factor must be a scalar/vector of positive integer(s)')
end

% number of factors
nf=length(factor);

% overall factor
of=prod(factor);

% check spacing
if(any(~strcmp(glgc(data,'leven'),'true')))
    error('seizmo:stretch:evenlySpacedOnly',...
        'illegal operation on unevenly spaced data');
end

% header info
[delta]=gh(data,'delta');
delta=delta*of;

% decimate and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % loop over components
    [len,ncmp]=size(data(i).x);
    newlen=(len-1)*of+1;
    save=zeros(newlen,ncmp);
    for j=1:ncmp
        temp=data(i).x(:,j);
        % loop over factors
        for k=1:nf
            % stretch (filter length 4, cutoff freq is nyquist)
            temp=interp(temp,factor(k),4,1);
        end
        
        % Matlab extrapolates past the end, so here
        % we truncate the extrapolated values off
        save(:,j)=temp(1:newlen,1);
    end
    
    % change class back
    data(i).x=oclass(save);
    
    % update header
    data(i)=ch(data(i),'npts',newlen,'delta',delta(i),...
        'depmin',min(data(i).x(:)),...
        'depmen',mean(data(i).x(:)),...
        'depmax',max(data(i).x(:)));
end

end
