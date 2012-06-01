function [data,idx]=splitrecords(data)
%SPLITRECORDS    Split up components into separate records
%
%    Usage:    data=splitrecords(data)
%              [data,idx]=splitrecords(data)
%
%    Description:
%     DATA=SPLITRECORDS(DATA) returns a dataset with all multiple component
%     records (spectral record components included) separated into single
%     component records.  Header info is replicated to all new records.
%     All spectral records are changed to xy records.
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
%     % JOINRECORDS does not undo SPLITRECORDS! This returns only 1 record
%     % if it succeeds because it attempts to join everything together:
%     data=joinrecords(splitrecords(data));
%
%     % Undo SPLITRECORDS (even iftype change):
%     iftype=getheader(data,'iftype');
%     [splitdata,idx]=splitrecords(data);
%     savedata=data;
%     for i=1:max(idx)
%         data(i)=joinrecords(splitdata(idx==i));
%     end
%     data=changeheader(data,'iftype',iftype);
%     % Now subtract and plot to make sure:
%     plot1(subtractrecords(savedata,data));
%
%    See also: JOINRECORDS, MULTIFUN, GETSPECTRALCMP, KEEPAM, KEEPPH,
%              KEEPRL, KEEPIM, CUT

%     Version History:
%        June 29, 2009 - initial version
%        Jan. 30, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb.  3, 2010 - fixed verbose msg bug, fixed bug when no mcmp
%        Jan.  6, 2011 - recordfun/multifun rename
%        Nov.  2, 2011 - doc update
%        Mar. 13, 2012 - use getheader improvements, fix example
%        May  31, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2012 at 16:40 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt splitting records
try
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    idx=(1:nrecs).';

    % get filetype and ncmp
    [ncmp,iftype]=getheader(data,'ncmp','iftype id');

    % double ncmp for spectral
    isspectral=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    if(any(isspectral)); ncmp(isspectral)=ncmp(isspectral)*2; end

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
    
    % detail message
    if(verbose)
        disp('Splitting Multi-Component Record(s)');
        print_time_left(0,nrecs);
    end

    % split records
    depmen=nan(totrecs,1); depmin=depmen; depmax=depmen;
    count=0;
    for i=mcmp
        for j=1:ncmp(i)
            count=count+1;
            idx2(count)=i;
            if(isspectral(i))
                iftype2(count)={'ixy'};
            else
                iftype2(count)=iftype(i);
            end
            for k=fieldsnodep.'
                data2(count).(k{:})=data(i).(k{:});
            end
            data2(count).dep=data(i).dep(:,j);
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % detail message
    if(verbose && (isempty(i) || i~=nrecs))
        print_time_left(nrecs,nrecs);
    end 

    % update header
    if(totrecs)
        data2=changeheader(data2,'ncmp',1,'iftype',iftype2,...
            'depmen',depmen,'depmin',depmin,'depmax',depmax);
    end
    
    % crunch/combine datasets
    data(mcmp)=[];
    data=[data(:); data2(:)];
    idx(mcmp)=[];
    idx=[idx; idx2];

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
