function [data]=correlate(data1,varargin)
%CORRELATE    Compute cross correlograms of SEIZMO data records
%
%    Description: CORRELATE(DATA) cross correlates every possible pairing
%     of records in DATA.  Output is a SEIZMO dataset of the correlograms.
%     
%     CORRELATE(DATA1,DATA2) cross correlates every record in DATA1 against
%     every record in DATA2 (meaning that records in DATA1 are the 'master'
%     records).
%
%     CORRELATE(...,'normxc',false) leaves the correlograms unnormalized
%     (they are normalized by their autocorrelations by default).
%
%    Header Changes:  B, E, DEPMEN, DEPMAX, DEPMIN, KT0-3
%     The master record for each correlation is used for remaining header
%     values.  Slave record's KNETWK, KSTNM, KHOLE, KCMPNM are pushed into
%     fields KT0-3.  B, E are the lag range in seconds adjusted for the
%     differential beginning times of the records. All records are checked
%     that they share the same reference time and delta.
%
%    Usage: correlograms=correlate(data)
%           correlograms=correlate(data,'normxc',false)
%           correlograms=correlate(data1,data2)
%           correlograms=correlate(data1,data2,'normxc',false)
%
%    Examples:
%     Equivalent to 'correlate' in SAC:
%     correlograms=correlate(data(1),data,'normxc',false)
%
%    See also: convolve, dft, idft

% split based on npeaks (0 or 1+)           check
% check datsets                             check
% peaks get added to fields .cg .lg .pg     check
% correlograms replace dataset
% figure out lags
%   lag_slave + b_slave - b_master
% - add fields to know
%   - master record
%   - slave record
% - change name to CORR_MASTER_SLAVE

% check nargin
if(nargin<1)
    error('seizmo:selectrecords:notEnoughInputs',...
        'Not enough input arguments.');
end

% check data structure
msg=seizmocheck(data1,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% check headers
data1=checkheader(data1);

% 1 or 2 datasets
nrecs1=numel(data1); onedata=true;
if(nargin>1 && isseizmo(varargin{1},'dep'))
    disp('yikes')
    data2=varargin{1};
    varargin(1)=[];
    nrecs2=numel(data2); 
    onedata=false;
    data2=checkheader(data2);
end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% check varargin
nvarargin=numel(varargin);
if(mod(nvarargin,2))
    error('seizmo:correlate:OptionMustBePaired',...
        'Options must be paired with a value!');
end

% how many peaks
npeaks=0; % default is no peak picking
for i=1:2:nvarargin
    if(ischar(varargin{i}) && strcmpi('npeaks',varargin{i}))
        npeaks=varargin{i+1};
    end
end

% check npeaks
if(~isscalar(npeaks) && ~isnumeric(npeaks) && fix(npeaks)~=npeaks)
    error('seizmo:correlate:badInput',...
        'Option NPEAKS must be a scalar integer')
end

% split based on npeaks
% npeaks>0 ==> assign grids to new fields
if(npeaks)
    % split based on ndatasets
    if(onedata)
        % check dataset
        delta=cross_check_data(data1);
        % get relative start times (disregarding absolute timing)
        b1=getheader(data1,'b');
        bdiff=(b1(:,ones(nrecs1,1)).'-b1(:,ones(nrecs1,1))).';
        bdiff=bdiff(tril(true(nrecs1),-1)); % extract upper triangle
        % extract records
        data1=combo(data1);
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(ndims(s):4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    else % two datasets
        % check datasets
        delta=cross_check_data(data1,data2);
        % get relative start times (disregarding absolute timing)
        b1=gh(data1,'b');
        b2=gh(varargin{1},'b');
        bdiff=b2(:,ones(nrecs1,1)).'-b1(:,ones(nrecs2,1));
        % extract records
        data1=combo(data1);
        varargin{1}=combo(varargin{1});
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(ndims(s):4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    end
% npeaks==0 ==> replace dataset with correlograms
else
    % split based on ndatasets
    if(onedata)
        % check datasets
        delta=cross_check_data(data1);
        
        % extract header info needed
        [knetwk,kstnm,khole,kcmpnm,b,stla,stlo,stdp,stel]=...
            getheader(data1,'knetwk','kstnm','khole','kcmpnm','b',...
            'stla','stlo','stdp','stel');
        
        % extract records
        data1=combinerecords(data1);
        
        % get correlation peaks
        [cg,lg]=mcxc(data1,varargin{:});
        
        % separate records
        [ncors,one,nlags]=size(cg);
        if(ncors==1)
            cg=squeeze(cg);
        else
            cg=squeeze(cg).';
        end
        cg=mat2cell(cg,nlags,ones(ncors,1));
        
        % make lags for each record
        % 1 to nrecs-1 outer loop is master
        % 2 to nrecs inner loop is slave
        bdiff=b(:,ones(nrecs1,1)).'-b(:,ones(nrecs1,1));
        bdiff=bdiff(tril(true(nrecs1),-1));
        lg=lg(:,ones(ncors,1))*delta-bdiff(:,ones(nlags,1)).';
        lg=mat2cell(lg,nlags,ones(ncors,1));
        
        % make a new dataset
        data=[lg; cg];
        data=bseizmo(data{:});
        
        % get names
        i=1:nrecs1;
        idx=cellstr(num2str(i.',...
            ['%0' num2str(ceil(log10(nrecs1+1))) 'd']));
        idx=idx(:,ones(nrecs1,1)).';
        knetwk=knetwk(:,ones(nrecs1,1));
        kstnm=kstnm(:,ones(nrecs1,1));
        khole=khole(:,ones(nrecs1,1));
        kcmpnm=kcmpnm(:,ones(nrecs1,1));
        names=strcat({'CORR_-_MASTER_-_REC'},idx,{'_-_'},...
            knetwk',{'.'},kstnm',{'.'},khole',{'.'},kcmpnm',...
            {'_-_SLAVE_-_REC'},idx',{'_-_'},...
            knetwk,{'.'},kstnm,{'.'},khole,{'.'},kcmpnm);
        names=names(tril(true(nrecs1),-1));
        [data.name]=deal(names{:});
        
        % add record numbers to header
        midx=i(ones(nrecs1,1),:);
        sidx=midx.';
        midx=midx(tril(true(nrecs1),-1));
        sidx=sidx(tril(true(nrecs1),-1));
        
        % add station locations to header
        evla=stla(:,ones(nrecs1,1)); stla=evla.';
        evlo=stlo(:,ones(nrecs1,1)); stlo=evlo.';
        evel=stel(:,ones(nrecs1,1)); stel=evel.';
        evdp=stdp(:,ones(nrecs1,1)); stdp=evdp.';
        stla=stla(tril(true(nrecs1),-1)); evla=evla(tril(true(nrecs1),-1));
        stlo=stlo(tril(true(nrecs1),-1)); evlo=evlo(tril(true(nrecs1),-1));
        stel=stel(tril(true(nrecs1),-1)); evel=evel(tril(true(nrecs1),-1));
        stdp=stdp(tril(true(nrecs1),-1)); evdp=evdp(tril(true(nrecs1),-1));
        
        % update header
        data=changeheader(data,...
            'user0',midx,'kuser0','MASTER',...
            'user1',sidx,'kuser1','SLAVE',...
            'stla',stla,'stlo',stlo,'stel',stel,'stdp',stdp,...
            'evla',evla,'evlo',evlo,'evel',evel,'evdp',evdp);
        % still need to update dep* !!!
        % filetype is not sac!!!
    else
        % first set is master trace
        % second is slave
        % lag*delta-(master-slave)
    end
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end

function [delta]=cross_check_data(data1,data2)
% assure datasets are capable of being sensibly cross correlated

% if 2 datasets, datasets must be consistent with one another
if(nargin==2)
    delta=getheader([data1(:); data2(:)],'delta');
    iftype=getenumdesc([data1(:); data2(:)],'iftype');
    leven=getlgc([data1(:); data2(:)],'leven');
    ncmp=getncmp([data1(:); data2(:)]);
% if 1 dataset, needs to be consistent with itself
else
    [delta]=getheader(data1,'delta');
    iftype=getenumdesc(data1,'iftype');
    ncmp=getncmp(data1);
    leven=getlgc(data1,'leven');
end

% check delta,iftype,leven,ncmp
if(~isscalar(unique(delta)))
    error('seizmo:correlate:mixedSampleRates',...
        'Mixed sample rates not allowed!');
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:correlate:illegalOperation',...
        'Illegal operation on spectral/xyz record(s)!');
end
if(any(~strcmpi(leven,'true')))
    error('seizmo:correlate:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced data!');
end
if(any(ncmp>1))
    error('seizmo:correlate:tooManyComponents',...
        'Multi-component data is not supported!')
end

% return common sample rate
delta=delta(1);

end
