function [data]=squish(data,factor)
%SQUISH    Downsample SEIZMO records by an integer factor
%
%    Usage:    data=squish(data,factors)
%
%    Description: SQUISH(DATA,FACTOR) downsamples SEIZMO records by an 
%     integer FACTOR after implementing an anti-aliasing filter.  To avoid
%     adding significant numerical noise to the data, keep the decimation
%     factor below 13.  If a decimation factor greater than this is needed,
%     consider using a cascade of decimation factors (by giving an array of
%     factors for FACTOR).  For a usage case, see the examples below.
%     Strong amplitudes at or near the start and end of records will
%     introduce edge effects that can be avoided by first detrending and 
%     then tapering records.  Uses the Matlab function decimate (Signal 
%     Processing Toolbox).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA, NPTS, E, DEPMEN, DEPMIN, DEPMAX
%
%    Examples: 
%     To halve the samplerates of records in data:
%      data=squish(data,2)
%
%     To cascade records to a samplerate 40 times lower;
%      data=squish(data,[5 8])
%
%    See also: stretch, syncrates, interpolate

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
%        Nov. 23, 2008 - update for new name schema (now SQUISH)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:00 GMT

% todo:

% check inputs
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

% empty factor
if(isempty(factor)); return; end

% check factor
if(~isnumeric(factor) || any(factor<1) || any(fix(factor)~=factor))
    error('seizmo:squish:badInput',...
        'FACTOR must be one or more positive integers!')
end

% number of factors
factor=factor(:);
nf=numel(factor);

% overall factor
of=prod(factor);

% check filetype
iftype=getenumdesc(data,'iftype');
if(strcmpi(iftype,'General XYZ (3-D) file'))
    error('seizmo:squish:illegalFiletype',...
        'Illegal operation on xyz records!');
elseif(any(strcmpi(iftype,'Spectral File-Real/Imag') | ...
        strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:squish:illegalFiletype',...
        'Illegal operation on spectral records!');
end

% check spacing
if(any(~strcmpi(getlgc(data,'leven'),'true')))
    error('seizmo:squish:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced records!');
end

% header info
[b,delta]=getheader(data,'b','delta');
delta=delta*of;

% number of records
nrecs=numel(data);

% decimate and update header
e=nan(nrecs,1); npts=e; 
depmen=e; depmin=e; depmax=e;
for i=1:nrecs
    % dataless support
    [npts(i),ncmp]=size(data(i).dep);
    if(any([npts(i) ncmp]==0)); npts(i)=0; continue; end
    
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
            e(i)=b(i)+(npts(i)-1)*delta(i);
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
data(i)=changeheader(data(i),'npts',npts,'delta',delta,'e',e,...
    'depmin',depmin,'depmen',depmen,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
