function [data]=addf(varargin)
%ADDF    Add SAClab data records
%
%    Description: ADDF(DATA) will add all records in DATA and will return a
%     single record with its header fields set to those of the first record
%     in DATA.  The header can be set to that of the last record by setting
%     option 'newhdr' to TRUE.  Records should be of the same filetype, be
%     evenly sampled, have the same sample rate, number of points, and 
%     timing but these can all be ignored (for better or for worse) by 
%     setting options available in BINOPERR to 'ignore'.
%     
%     ADDF(DATA1,DATA2) will add the records in DATA2 to DATA1.  If either
%     DATA1 or DATA2 are a single record the record will be added to every
%     record of the other dataset.  DATA1 and DATA2 must contain the same
%     number of records otherwise.
%     
%     ADDF(DATA1,DATA2,...,DATAN) adds records in all N datasets such that
%     the first record of each dataset is added together and so on.
%     Therefore every dataset must have the same number of records.  The 
%     exception to this rule is single record datasets, which are
%     replicated to match the size of the rest of the datasets.
%     
%     ADDF(...,'newhdr',true) controls the inheritance of header fields and
%     will set the resultant records' header fields to those of the last 
%     record to be added on.  So 
%                    ADDF(DATA1,DATA2,'newhdr',true)
%     will produce records with header fields set to those in DATA2.  By
%     default 'newhdr' is set to FALSE which sets the resultant records'
%     header fields to those in DATA1.  If adding all records in a single
%     dataset, setting 'newhdr' to TRUE will set the resultant record's 
%     header equal to the last record's header.  Leaving 'newhdr' set to
%     the default FALSE will set the resultant record's header to that of
%     the first record's header.
%
%     The following options may also be controlled by BINOPERR.
%     
%     ADDF(...,'npts','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records with different numbers of points.  If the
%     option is set to 'warn' or 'ignore', the number of points in the
%     resultant records will be equal to that of the shortest record.  
%     Note that points are added according to their order in the record not
%     by their timing, such that the first points are always added together
%     and so on.  By default 'npts' is set to 'error'.
%     
%     ADDF(...,'delta','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records with different sample rates.  If the option is
%     set to 'warn' or 'ignore', the records are just added point for
%     point (basically ignoring timing).  The resultant records' sample 
%     rates are determined by the parent of their header fields (set by 
%     option 'newhdr').  By default 'delta' is set to 'error'.
%     
%     ADDF(...,'begin','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records with different begin times.  If the option is
%     set to 'warn' or 'ignore', the resultant records' begin times are 
%     determined by the parent of their header fields (set by option 
%     'newhdr').  By default 'begin' is set to 'warn'.
%     
%     ADDF(...,'ref','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records with different reference times.  If the option
%     is set to 'warn' or 'ignore', the resultant records' reference times
%     are determined by the parent of their header fields (set by option 
%     'newhdr').  By default 'ref' is set to 'warn'.
%     
%     ADDF(...,'ncmp','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records with different numbers of components.  If the
%     option is set to 'warn' or 'ignore', the number of components in the
%     resultant records will be equal to that of the record with the least.  
%     Note that components are added according to their order in the record
%     so that the first components always add.  By default 'ncmp' is set to
%     'error'.
%     
%     ADDF(...,'leven','error|warn|ignore') sets the behavior of how ADDF
%     handles adding unevenly sampled records.  If the option is set to 
%     'warn' or 'ignore', the records are just added point for point 
%     (basically ignoring timing).  The resultant records' leven fields are
%     determined by the parent of their header fields (set by option 
%     'newhdr').  By default 'leven' is set to 'error'.
%     
%     ADDF(...,'iftype','error|warn|ignore') sets the behavior of how ADDF
%     handles adding records of different types.  If the option is set to 
%     'warn' or 'ignore', the records are just added point for point.  The 
%     resultant records' iftypes are determined by the parent of their 
%     header fields (set by option 'newhdr').  By default 'iftype' is set 
%     to 'error'.
%     
%    Notes:
%     - Ampl-Phase spectral records are first converted to Real-Imag to
%       assure the operation is linear and equal to that on Real-Imag 
%       records.  If you want to workaround this, convert the Ampl-Phase
%       records to General X vs Y.
%
%    System requirements: Matlab 7
%
%    Data requirements: Records to be added should match in NPTS, DELTA, B,
%     IFTYPE, number of components and reference time.  These can be 
%     ignored by setting options with BINOPERR or inline options.
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN,
%     NPTS, E, NCMP (see option 'npts' and 'ncmp')
%     See option 'newhdr' for inheritance of other header fields.
%    
%    Usage: [data]=addf(data)
%           [data]=addf(data1,data2)
%           [data]=addf(data1,data2,...,dataN)
%           [data]=addf(...,'newhdr',true)
%           [data]=addf(...,'npts','error|warn|ignore')
%           [data]=addf(...,'delta','error|warn|ignore')
%           [data]=addf(...,'begin','error|warn|ignore')
%           [data]=addf(...,'ref','error|warn|ignore')
%           [data]=addf(...,'ncmp','error|warn|ignore')
%           [data]=addf(...,'leven','error|warn|ignore')
%           [data]=addf(...,'iftype','error|warn|ignore')
%
%    Examples:
%     Display a stack of the records:
%      p1(addf(data))
%
%     Add records from one dataset to another
%      data=addf(data1,data2)
%
%    See also: subf, mulf, divf, binoperr

%     Version History:
%        June 10, 2008 - initial version
%        June 11, 2008 - full filetype and class support
%        June 20, 2008 - documentation update, 'ncmp' option
%        June 28, 2008 - fixed amph2rlim handling, documentation update,
%                        .dep and .ind rather than .x and .t
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 28, 2008 at 21:45 GMT

% todo:
% - dataless support

% default options
option.NEWHDR=false;
option.NPTS='error';
option.DELTA='error';
option.BEGIN='warn';
option.REF='warn';
option.NCMP='error';
option.LEVEN='error';
option.IFTYPE='error';

% get options set by BINOPERR (SACLAB global)
global SACLAB; fields=fieldnames(option).';
for i=fields; if(isfield(SACLAB,i)); option.(i{:})=SACLAB.(i{:}); end; end

% find all datasets in inline arguments
isdata=false(1,nargin);
for i=1:nargin; isdata(i)=isseis(varargin{i},'dep'); end

% push datasets into a separate variable
data=varargin(isdata);
varargin(isdata)=[];

% options must be field-value pairs
nargopt=length(varargin);
if(mod(nargopt,2))
    error('SAClab:addf:badNumOptions','Unpaired option!')
end

% get inline options
for i=1:2:nargopt
    % plotting parameter field
    varargin{i}=upper(varargin{i});
    if(isfield(option,varargin{i}))
        option.(varargin{i})=varargin{i+1};
    else
        warning('SAClab:addf:badInput','Unknown Option: %s',varargin{i}); 
    end
end

% get number of records in each dataset
ndatasets=length(data);
nrecs=zeros(1,ndatasets);
for i=1:ndatasets
    nrecs(i)=numel(data{i});
end

% check for bad sized datasets
maxrecs=max(nrecs);
if(any(nrecs~=1 & nrecs~=maxrecs))
    error('SAClab:addf:nrecsMismatch','Number of records inconsistent')
end

% expand scalar datasets
for i=find(nrecs==1)
    data{i}(1:maxrecs)=data{i};
end

% get header fields
b=cell(1,ndatasets); npts=b; delta=b; leven=b; iftype=b; ncmp=b;
nzyear=b; nzjday=b; nzhour=b; nzmin=b; nzsec=b; nzmsec=b;
for i=1:ndatasets
    leven{i}=glgc(data{i},'leven');
    iftype{i}=genumdesc(data{i},'iftype');
    [npts{i},delta{i},b{i},nzyear{i},nzjday{i},...
        nzhour{i},nzmin{i},nzsec{i},nzmsec{i}]=...
        gh(data{i},'npts','delta','b',...
        'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec');
    ncmp{i}=zeros(maxrecs,1);
    for j=1:maxrecs; ncmp{i}(j)=size(data{i}(j).dep,2); end
end

% 2+ datasets
if(ndatasets>1)
    % check records
    if(~isequal(iftype{:}))
        report.identifier='SAClab:addf:mixedIFTYPE';
        report.message='Filetypes differ for some records';
        if(strcmpi(option.IFTYPE,'error')); error(report);
        elseif(strcmpi(option.IFTYPE,'warn')); warning(report.identifier,report.message);
        end
    end
    for i=1:ndatasets
        if(any(~strcmpi(leven{i},'true')))
            report.identifier='SAClab:addf:illegalOperation';
            report.message='illegal operation on unevenly spaced record';
            if(strcmpi(option.LEVEN,'error')); error(report);
            elseif(strcmpi(option.LEVEN,'warn')); warning(report.identifier,report.message);
            end
        end
    end
    if(~isequal(ncmp{:}))
        report.identifier='SAClab:addf:mixedNCMP';
        report.message='Number of components differ for some records';
        if(strcmpi(option.NCMP,'error')); error(report);
        elseif(strcmpi(option.NCMP,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isequal(npts{:}))
        report.identifier='SAClab:addf:mixedNPTS';
        report.message='Number of points differ for some records';
        if(strcmpi(option.NPTS,'error')); error(report);
        elseif(strcmpi(option.NPTS,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isequal(delta{:}))
        report.identifier='SAClab:addf:mixedDELTA';
        report.message='Sample rates differ for some records';
        if(strcmpi(option.DELTA,'error')); error(report);
        elseif(strcmpi(option.DELTA,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isequal(b{:}))
        report.identifier='SAClab:addf:mixedB';
        report.message='Begin times differ for some records';
        if(strcmpi(option.BEGIN,'error')); error(report);
        elseif(strcmpi(option.BEGIN,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isequal(nzyear{:}) || ~isequal(nzjday{:}) ...
            || ~isequal(nzhour{:}) || ~isequal(nzmin{:}) ...
            || ~isequal(nzsec{:})  || ~isequal(nzmsec{:}))
        report.identifier='SAClab:addf:mixedReferenceTimes';
        report.message='Reference times differ for some records';
        if(strcmpi(option.REF,'error')); error(report);
        elseif(strcmpi(option.REF,'warn')); warning(report.identifier,report.message);
        end
    end
    
    % save class and convert to double precision
    oclass=cell(1,maxrecs);
    if(option.NEWHDR)
        for i=1:maxrecs; oclass{i}=str2func(class(data{end}(i).dep)); end
    else
        for i=1:maxrecs; oclass{i}=str2func(class(data{1}(i).dep)); end
    end
    for i=1:ndatasets
        for j=1:maxrecs; data{i}(j).dep=double(data{i}(j).dep); end
    end
    
    % convert amplitude-phase files to real-imaginary so that operations
    % are consistent. amph2rlim must be done here (before newhdr swap)
    if(option.NEWHDR)
        convertback=strcmpi(iftype{end},'Spectral File-Ampl/Phase');
    else
        convertback=strcmpi(iftype{1},'Spectral File-Ampl/Phase');
    end
    for i=1:ndatasets
        convert=strcmpi(iftype{i},'Spectral File-Ampl/Phase');
        if(any(convert)); data{i}(convert)=amph2rlim(data{i}(convert)); end
    end
    
    % newhdr flag (swap first and last dataset)
    if(option.NEWHDR); data([end 1])=data([1 end]); end
    
    % get min npts and ncmp in each added set
    allnpts=cell2mat(npts);
    minpts=min(allnpts,[],2);
    allncmp=cell2mat(ncmp);
    mincmp=min(allncmp,[],2);
    
    % add records
    for i=1:maxrecs
        for j=2:ndatasets
            data{1}(i).dep=...
                data{1}(i).dep(1:minpts(i),1:mincmp(i))...
                +data{j}(i).dep(1:minpts(i),1:mincmp(i));
        end
        
        % trim t field for unevenly spaced files
        if(isfield(data{1}(i),'ind') && ~isempty(data{1}(i).ind))
            data{1}(i).ind=data{1}(i).ind(1:minpts);
        end
        
        % change class back
        data{1}(i).dep=oclass{i}(data{1}(i).dep);
    end
    
    % update header
    data=chkhdr(data{1});
    
    % convert back to amph if necessary
    if(any(convertback))
        data(convertback)=rlim2amph(data(convertback));
    end
% 1 dataset
else
    % uncell
    data=data{:};
    
    % no records to add on
    if(isscalar(data)); return; end
    
    % check records
    if(~isscalar(unique(iftype{:})))
        report.identifier='SAClab:addf:mixedIFTYPE';
        report.message='Filetypes differ for some records';
        if(strcmpi(option.IFTYPE,'error')); error(report);
        elseif(strcmpi(option.IFTYPE,'warn')); warning(report.identifier,report.message);
        end
    end
    if(any(~strcmpi(leven{:},'true')))
        report.identifier='SAClab:addf:illegalOperation';
        report.message='illegal operation on unevenly spaced record';
        if(strcmpi(option.LEVEN,'error')); error(report);
        elseif(strcmpi(option.LEVEN,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isscalar(unique(ncmp{:})))
        report.identifier='SAClab:addf:mixedNCMP';
        report.message='Number of components differ for some records';
        if(strcmpi(option.NCMP,'error')); error(report);
        elseif(strcmpi(option.NCMP,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isscalar(unique(npts{:})))
        report.identifier='SAClab:addf:mixedNPTS';
        report.message='Number of points differ for some records';
        if(strcmpi(option.NPTS,'error')); error(report);
        elseif(strcmpi(option.NPTS,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isscalar(unique(delta{:})))
        report.identifier='SAClab:addf:mixedDELTA';
        report.message='Sample rates differ for some records';
        if(strcmpi(option.DELTA,'error')); error(report);
        elseif(strcmpi(option.DELTA,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isscalar(unique(b{:})))
        report.identifier='SAClab:addf:mixedB';
        report.message='Begin times differ for some records';
        if(strcmpi(option.BEGIN,'error')); error(report);
        elseif(strcmpi(option.BEGIN,'warn')); warning(report.identifier,report.message);
        end
    end
    if(~isscalar(unique(nzyear{:})) || ~isscalar(unique(nzjday{:})) || ...
            ~isscalar(unique(nzhour{:})) || ~isscalar(unique(nzmin{:})) ...
            || ~isscalar(unique(nzsec{:})) || ~isscalar(unique(nzmsec{:})))
        report.identifier='SAClab:addf:mixedReferenceTimes';
        report.message='Reference times differ for some records';
        if(strcmpi(option.REF,'error')); error(report);
        elseif(strcmpi(option.REF,'warn')); warning(report.identifier,report.message);
        end
    end
    
    % save class and convert to double precision
    if(option.NEWHDR); oclass=str2func(class(data(end).dep)); 
    else oclass=str2func(class(data(1).dep));
    end
    for i=1:nrecs; data(i).dep=double(data(i).dep); end
    
    % convert amplitude-phase files to real-imaginary so that operations
    % are consistent. amph2rlim must be done here (before newhdr)
    if(option.NEWHDR)
        convertback=strcmpi(iftype{:}(end),'Spectral File-Ampl/Phase');
    else
        convertback=strcmpi(iftype{:}(1),'Spectral File-Ampl/Phase');
    end
    convert=strcmpi(iftype{:},'Spectral File-Ampl/Phase');
    if(any(convert)); data(convert)=amph2rlim(data(convert)); end
    
    % newhdr flag (swap first and last record)
    if(option.NEWHDR); data([end 1])=data([1 end]); end
    
    % add records
    minpts=min(npts{:});
    mincmp=min(ncmp{:});
    for i=2:nrecs
        data(1).dep=data(1).dep(1:minpts,1:mincmp)...
                 +data(i).dep(1:minpts,1:mincmp);
    end
    
    % reduce to first record
    data=data(1);
    
    % trim ind field for unevenly spaced files
    if(isfield(data,'ind') && ~isempty(data.ind))
        data.ind=data.ind(1:minpts);
    end
    
    % change class back
    data.dep=oclass(data.dep);
    
    % update header
    data=chkhdr(data);
    
    % convert back to amph if necessary
    if(convertback); data=rlim2amph(data); end
end

end
