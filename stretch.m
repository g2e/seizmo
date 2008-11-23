function [data]=stretch(data,factor)
%STRETCH    Upsample SEIZMO records by an integer factor
%
%    Description: STRETCH(DATA,FACTOR) upsamples SEIZMO records by an 
%     integer FACTOR.  Anti-aliasing is not an issue so there is no limit 
%     imposed on the stretch factor.  Cascades are allowed regardless.  
%     Uses the Matlab function interp (Signal Processing Toolbox).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=stretch(data,factor)
%
%    Examples: 
%     To double samplerates:
%      data=stretch(data,2)
%
%     To cascade to a samplerate 40 times higher:
%      data=stretch(data,[5 8])
%
%    See also: squish, syncrates, interpolate

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 30, 2008 - doc update, code cleaning
%        Mar.  4, 2008 - major code cleaning, fix extrapolation
%        May  12, 2008 - fix dep* formula
%        June 15, 2008 - doc update, history added
%        Nov. 23, 2008 - doc update, history fixed, better checking, single
%                        changeheader call, update for new name schema,
%                        .dep rather than .x
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 23, 2008 at 09:45 GMT

% todo:

% check inputs
error(nargchk(2,2,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% empty factor
if(isempty(factor)); return; end

% check factor
if(~isnumeric(factor) || any(factor<1) || any(fix(factor)~=factor))
    error('seizmo:stretch:badInput',...
        'FACTOR must be one or more positive integers!');
end

% number of factors
factor=factor(:);
nf=numel(factor);

% overall factor
of=prod(factor);

% check filetype
iftype=getenumdesc(data,'iftype');
if(strcmpi(iftype,'General XYZ (3-D) file'))
    error('seizmo:deci:illegalFiletype',...
        'Illegal operation on xyz records!');
elseif(any(strcmpi(iftype,'Spectral File-Real/Imag') | ...
        strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:deci:illegalFiletype',...
        'Illegal operation on spectral records!');
end

% check spacing
if(any(~strcmpi(getlgc(data,'leven'),'true')))
    error('seizmo:stretch:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced records!');
end

% header info
delta=getheader(data,'delta');
delta=delta*of;

% number of records
nrecs=numel(data);

% decimate and update header
npts=nan(nrecs,1); depmen=npts; depmin=npts; depmax=npts;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); npts(i)=0; continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).x=double(data(i).dep);
    
    % loop over components (b/c interp can't handle arrays)
    [len,ncmp]=size(data(i).dep);
    npts(i)=(len-1)*of+1;
    save=zeros(npts(i),ncmp);
    for j=1:ncmp
        temp=data(i).dep(:,j);
        % loop over factors
        for k=1:nf
            % stretch (filter length 4, cutoff freq is nyquist)
            temp=interp(temp,factor(k),4,1);
        end
        
        % Matlab extrapolates past the end, so here
        % we truncate the extrapolated values off
        save(:,j)=temp(1:npts(i),1);
    end
    
    % change class back
    data(i).dep=oclass(save);
    
    % get dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=changeheader(data,'npts',npts,'delta',delta,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
