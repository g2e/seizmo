function [bytes]=seissize(data)
%SEISSIZE    Returns size of seismic file in bytes
%
%    Description: Calculates the expected size of a binary seismic
%     datafile in bytes using the input seislab structure.
%
%    Usage:  bytes=seissize(data)
%
%    See also: rdata

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% pull necessary header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
warning('off','seislab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','seislab:gh:fieldInvalid')

% check leven
t=strcmp(leven,'true');
f=strcmp(leven,'false');
if(~all(t | f))
    error('sieslab:seissize:levenBad',...
        'logical field leven needs to be set'); 
end

% fix and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('seislab:seissize:badNumCmp',...
        'field ncmp must be a positive integer')
end

% headers setup
v=[data.version].';
vers=unique(v);
nver=length(vers);
h(nver)=seishi(vers(nver));
for i=1:nver-1
    h(i)=seishi(vers(i));
end
nrecs=length(data);
v2=sum((v(:,ones(1,nver))==vers(:,ones(1,nrecs)).').*...
    repmat(1:nver,nrecs,1),2);
hdata=[h(v2).data];

% filetype
count=strcmp(iftype,'Time Series File')...
    +strcmp(iftype,'General X vs Y file')...
    +strcmp(iftype,'General XYZ (3-D) file')...
    +2*strcmp(iftype,'Spectral File-Real/Imag')...
    +2*strcmp(iftype,'Spectral File-Ampl/Phase');

% multi-component
count=ncmp.*count;

% uneven sampling
count=count+f;

% final tally
bytes=[hdata.startbyte].'+count.*npts.*[hdata.bytesize].';

end
