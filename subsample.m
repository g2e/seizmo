function [data]=deci(data,factor)
%DECI    Downsample SEIZMO data records by integer factor
%
%    Description: DECI(DATA,FACTOR) downsamples SEIZMO data records by an 
%     integer FACTOR after implementing an anti-aliasing filter.  To avoid
%     adding significant numerical noise to the data, keep the decimation
%     factor below 13.  If a decimation factor greater than this is needed,
%     consider using a cascade of decimation factors (by giving a vector
%     for FACTOR).  For a usage case for this see the examples below.
%     Strong amplitudes at or near the start and end of records will
%     introduce edge effects that can be reduced by first detrending and 
%     then tapering records.  Uses the Matlab function decimate (Signal 
%     Processing Toolbox).
%
%    Notes:
%
%    System requirements: Matlab 7, Signal Processing Toolbox
%
%    Data requirements: Evenly Spaced; Time Series or General X vs Y
%
%    Header changes: DELTA, NPTS, E, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:  data=deci(data,factors)
%
%    Examples: 
%     To halve the samplerates of records in data:
%      data=deci(data,2)
%
%     To cascade records to a samplerate 40 times lower;
%      data=deci(data,[5 8])
%
%    See also: stretch, syncsr, interpolate

%     Version History:
%        Oct. 30, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 30, 2008 - better input checking and doc update
%        Feb. 23, 2008 - fixed bugs on multi-component files, improved
%                        input checks, added examples
%        Feb. 28, 2008 - seischk support
%        May  12, 2008 - dep* fix
%        June 15, 2008 - doc update
%        June 30, 2008 - single ch call, dataless support, handle 
%                        decimation to single point, doc update,
%                        .dep rather than .x, filetype checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 30, 2008 at 06:20 GMT

% todo:
%

% check inputs
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'dep'))

% empty factor
if(isempty(factor)); return; end

% check factor
if(~isnumeric(factor) || ~isvector(factor) || ...
        any(factor<1) || any(fix(factor)~=factor))
    error('seizmo:deci:badInput',...
        'factor must be a scalar/vector of positive integer(s)')
end

% number of factors
nf=length(factor);

% overall factor
of=prod(factor);

% check filetype
iftype=genumdesc(data,'iftype');
if(strcmpi(iftype,'General XYZ (3-D) file'))
    error('seizmo:deci:illegalFiletype',...
        'illegal operation on xyz file');
elseif(any(strcmpi(iftype,'Spectral File-Real/Imag') | ...
        strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:deci:illegalFiletype',...
        'illegal operation on spectral file');
end

% check spacing
if(any(~strcmpi(glgc(data,'leven'),'true')))
    error('seizmo:deci:evenlySpacedOnly',...
        'illegal operation on unevenly spaced data');
end

% header info
[b,delta]=gh(data,'b','delta');
delta=delta*of;

% number of records
nrecs=numel(data);

% decimate and update header
e=nan(nrecs,1); npts=e; 
depmen=e; depmin=e; depmax=e;
for i=1:nrecs
    % dataless support
    [npts(i),ncmp]=size(data(i).dep);
    if(any([npts(i) ncmp]==0)); npts(i)=0; delta(i)=nan; continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % loop over components (because decimate can't handle arrays)
    for j=1:ncmp
        temp=data(i).dep(:,j);
        % loop over factors
        for k=1:nf
            temp=decimate(temp,factor(k));
        end
        % preallocate after we know length
        if(j==1)
            npts(i)=size(temp,1);
            if(npts(i)==1)
                delta(i)=nan;
                e(i)=b(i);
            else
                e(i)=b(i)+(npts(i)-1)*delta(i);
            end
            save=zeros(npts(i),ncmp);
        end
        save(:,j)=temp;
    end
    
    % change class back
    data(i).dep=oclass(save);
    
    % dep* stats
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update header
data(i)=ch(data(i),'npts',npts,'delta',delta,'e',e,'depmin',depmin,...
    'depmen',depmen,'depmax',depmax);

end
