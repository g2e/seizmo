function [data]=chkhdr(data,varargin)
%CHKHDR    Check and update header field/values of SAClab records
%
%    Description: CHKHDR(DATA) currently just updates the header fields 
%     'depmin', 'depmax', 'depmen', 'npts' and 'e' to reflect the records 
%     in DATA.  Field 'b' is also updated for unevenly sampled records.  
%     
%     CHKHDR(DATA,FIELD1,VALUE1,...,FIELDN,VALUEN) allows additional header
%     changes to be passed through to CH.
%
%    Notes:
%     - Checks that the LEVEN field is set
%
%    Requirements: Matlab 7
%
%    Supported types: ALL
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, E, B
%
%    Usage: data=chkhdr(data)
%           data=chkhdr(data,'field1',values1,...,'fieldN',valuesN)
%
%    Examples:
%     Check some basic header fields and also modify the t0 field with 
%     some arrival data:
%      data=chkhdr(data,'t0',P_arrivaltimes)
%
%    See also:  fixdelta, ch, gh

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Updated documentation
%        June 20, 2008 - Updates header then checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 20, 2008 at 07:30 GMT

% todo
%
% dataless check and update
% leven check and fix
% lat lon check and fix
% lcalda check and update
% spectral file exclusion for dep*
% 

% check input
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% change header first
data=ch(data,varargin{:});

% lat lon check
%[evla,evlo,stla,stlo]=gh(data,'evla','evlo','stla','stlo');

% lcalda check
%[lcalda]=glgc(data,'lcalda');

% leven check/fix
%leven=glgc(data,'leven');

% get header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
tru=strcmp(leven,'true');
fals=strcmp(leven,'false');

% number of records
nrecs=length(data);

% preallocate
len=zeros(nrecs,1); 
depmax=len; depmin=len; depmen=len;

% loop over all records, get npts, dep(min,max,men)
for i=1:nrecs
    len(i)=size(data(i).x,1);
    depmax(i)=max(data(i).x(:));
    depmin(i)=min(data(i).x(:));
    depmen(i)=mean(data(i).x(:));
end

% work only on uneven records, updating b and e fields
for i=find(fals).'
    tlen=size(data(i).t,1);
    if(tlen~=len(i))
        error('SAClab:chkhdr:dataSizeInconsistent',...
            ['number of samples in timing vector inconsistent with'...
            'data size for record %d'],i);
    end
    data(i)=ch(data(i),'b',data(i).t(1),'e',data(i).t(end));
end 

% work only on even records, updating e field
[b,delta]=gh(data(tru),'b','delta'); 
data(tru)=ch(data(tru),'e',b+(len(tru)-1).*delta);

% update some header fields for all records
data=ch(data,'depmax',depmax,'npts',len,...
    'depmin',depmin,'depmen',depmen);

end
