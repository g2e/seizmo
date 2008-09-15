function [bytes]=seissize(data)
%SEISSIZE    Returns header estimated size of SAClab datafiles in bytes
%
%    Description: SEISSIZE(DATA) returns the expected size (in bytes) using
%     the header info of each datafile that would be written from DATA. 
%     This is mainly used to rapidly detect files on disk that are 
%     inconsistent in size from what their header info indicates.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: expects SAClab DATA structure to have
%     version and head fields filled
%
%    Header changes: N/A
%
%    Usage:  bytes=seissize(data)
%
%    Examples:
%     Calculate expected filesizes for SAC files in current directory:
%      seissize(rh('*.SAC'))
%
%    See also: rh, rdata

%     Version History:
%        Feb. 29, 2008 - initial version (from subfunction sacsize)
%        Mar.  4, 2008 - doc update, use LGCCHK
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 21:15 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% pull necessary header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
warning('off','SAClab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% fix and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:seissize:badNumCmp',...
        'field ncmp must be a positive integer')
end

% headers setup
v=[data.version].';
vers=unique(v);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
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
count=count+strcmp(leven,'false');

% final tally
bytes=[hdata.startbyte].'+count.*npts.*[hdata.bytesize].';

end
