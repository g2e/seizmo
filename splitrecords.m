function [data,idx]=splitrecords(data)
%SPLITRECORDS    Split up components into separate records
%
%    Usage:    data=splitrecords(data)
%              [data,idx]=splitrecords(data)
%
%    Description: SPLITRECORDS(DATA) returns a dataset with all multiple
%     component records (spectral record components included) separated
%     into single component records.  Header info is replicated to all new
%     records.  All spectral records are changed to xy records.
%
%     [DATA,IDX]=SPLITRECORDS(DATA) also returns an index array to indicate
%     which records in the input dataset the output dataset records came
%     from.  So IDX will have the same number of elements as the output
%     dataset and MAX(IDX) is equal to the number of records in the input
%     dataset.  This is useful for undoing SPLITRECORDS.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NCMP, IFTYPE
%
%    Examples:
%     This returns only 1 record if it succeeds:
%      data=joinrecords(splitrecords(data));
%
%     Undo splitrecords (even iftype change):
%      [splitdata,idx]=splitrecords(data);
%      iftype=getheader(data,'iftype');
%      savedata=data;
%      for i=1:max(idx)
%          data(i)=joinrecords(splitdata(idx==i));
%      end
%      data=changeheader(data,'iftype',iftype);
%     Now subtract and plot to make sure:
%      p1(subtractrecords(savedata,data));
%
%    See also: joinrecords, recordfun, getspectralcmp, keepam, keepph,
%              keeprl, keepim, cut

%     Version History:
%        June 29, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 29, 2009 at 00:45 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% number of records
nrecs=numel(data);
idx=(1:nrecs).';

% get filetype and ncmp
iftype=getenumdesc(data,'iftype');
ncmp=getncmp(data);

% double ncmp for spectral
isspectral=strcmpi(iftype,'Spectral File-Real/Imag')...
    | strcmpi(iftype,'Spectral File-Ampl/Phase');
if(any(isspectral))
    ncmp(isspectral)=ncmp(isspectral)*2;
end

% total new records
mcmp=find(ncmp.'>1);
totrecs=sum(ncmp(mcmp));

% allocate
idx2=nan(totrecs,1);
iftype2=cell(totrecs,1);

% make a new struct
fields=fieldnames(data);
fieldsnodep=setxor(fields,'dep');
allocate=[fields.'; cell(1,numel(fields))];
data2(1:totrecs,1)=struct(allocate{:});

% split records
depmen=nan(totrecs,1); depmin=depmen; depmax=depmen;
count=0;
for i=mcmp
    for j=1:ncmp(i)
        count=count+1;
        idx2(count)=i;
        if(isspectral(i))
            iftype2(count)={'General X vs Y file'};
        else
            iftype2(count)=iftype(i);
        end
        for k=fieldsnodep.'
            data2(count).(k{:})=data(i).(k{:});
        end
        data2(count).dep=data(i).dep(:,j);
    end
end

% update header
warning('off','seizmo:changeheader:fieldInvalid')
data2=changeheader(data2,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'ncmp',1,'iftype',iftype2);
warning('on','seizmo:changeheader:fieldInvalid')

% crunch/combine datasets
data(mcmp)=[];
data=[data(:); data2(:)];
idx(mcmp)=[];
idx=[idx; idx2];

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
