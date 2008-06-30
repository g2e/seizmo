function [data]=dif(data)
%DIF    Differentiate SAClab data records using discrete differences
%
%    Description: DIF(DATA) calculates and returns the derivative of each 
%     record in the SAClab structure DATA using the differences between 
%     points as an approximation of the derivative at the midpoint.  Works 
%     with unevenly spaced data.
%
%    Notes: 
%     - Timing is shifted to midpoints
%     - Begin time increases by a half sample interval
%     - Reduces npts by one
%
%    System requirements: Matlab 7
%
%    Data requirements: Time Series or General X vs Y
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, B, E
%
%    Usage: data=dif(data)
%
%    Examples:
%     These are equal:
%      rtrend(data)
%      integrt2(dif(data))
%
%    See also: integrt, integrt2

%    Version History:
%       Jan. 28, 2008 - initial version
%       Feb. 23, 2008 - glgc support
%       Feb. 28, 2008 - seischk support
%       Mar.  4, 2008 - lgcchk support
%       May  12, 2008 - dep* fix
%       June 20, 2008 - minor documentation update
%       June 29, 2008 - documentation update, .dep & .ind rather than .x &
%                       .t, dataless support, only calls ch once, strict
%                       filetype check
%
%    Written by Garrett Euler (ggeuler at wustl dot edu)
%    Last Updated June 29, 2008 at 07:25 GMT

% todo:
% - 3pt & 5pt difference operator

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% check for unsupported filetypes
iftype=genumdesc(data,'iftype');
if(strcmpi(iftype,'General XYZ (3-D) file'))
    error('SAClab:dif:illegalFiletype',...
        'illegal operation on xyz file');
elseif(any(strcmpi(iftype,'Spectral File-Real/Imag') | ...
        strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:dif:illegalFiletype',...
        'illegal operation on spectral file');
end

% retreive header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
[delta,b,e,npts]=gh(data,'delta','b','e','npts');

% number of records
nrecs=numel(data);

% take derivative and update header
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % dataless support
    if(npts(i)==0); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % evenly spaced
    if(strcmpi(leven(i),'true'))
        data(i).dep=diff(data(i).dep)/delta(i);
        b(i)=b(i)+delta(i)/2; e(i)=e(i)-delta(i)/2; npts(i)=npts(i)-1;
    % unevenly spaced
    else
        data(i).ind=double(data(i).ind);
        ind=diff(data(i).ind);
        data(i).dep=diff(data(i).dep)./ind(:,ones(1,size(data(i).dep,2)));
        data(i).ind=oclass(data(i).ind(1:npts(i)-1)+ind/2);
        npts(i)=npts(i)-1; b(i)=data(i).ind(1); e(i)=data(i).ind(end);
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % get dep* info
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'b',b,'e',e,'npts',npts);

end
